# VM Sandboxes

<Callout variant="beta" />

Sandboxes can be run on top of a full virtual machine rather than on top of gVisor. This gives
each Sandbox a real Linux kernel, which makes certain workloads (e.g. Docker systems) behave
the way they would on a normal Linux host.

You can use the VM runtime for your Sandbox by passing `experimental_options={"vm_runtime": True}`
to `Sandbox.create()`.

<Collapsible title="VM demo">

```python fixture:sb_app
with modal.enable_output():
    sb = modal.Sandbox.create(
        app=sb_app,
        cpu=2,  # physical cores
        memory=4096,  # MiB
        experimental_options={"vm_runtime": True},
    )

# add a script that uses VM Sandbox features
sb.filesystem.write_text(
    """
  # Format an ext4 filesystem onto a regular file.
  truncate -s 100M /tmp/disk.img
  mkfs.ext4 -F /tmp/disk.img

  # Mount it. This works in a VM, but isn't supported in gVisor.
  mkdir -p /mnt/loop
  mount -o loop /tmp/disk.img /mnt/loop
""",
    "/tmp/mount_loopback_filesystem.sh",
)

p = sb.exec("bash", "/tmp/mount_loopback_filesystem.sh")
p.wait()

print(p.stdout.read())

print(p.stderr.read())

assert p.returncode == 0  # error if the program in the Sandbox fails

sb.terminate()
```

</Collapsible>

VM Sandboxes are also the recommended method to run Docker in Sandboxes. To try this out,
copy the following program to e.g. `docker_in_modal_demo.py`, and run it with
`python docker_in_modal_demo.py`.

<Collapsible title="Docker-in-Sandbox demo">

<!-- Keep the code block below in sync with synthetic_monitoring/benchmarks/docker_in_modal.py.
The "marker" comments below are used to diff this code with that in the synmon. Don't change them. -->

<!-- synmon-sync:docker_in_modal:begin -->

```python
import modal

# Create an image for the parent Modal Sandbox, with Docker installed.
def create_modal_sandbox_image():
    image = (
        modal.Image.from_registry("ubuntu:24.04")
        .env({"DEBIAN_FRONTEND": "noninteractive"})
        .apt_install(["docker.io", "docker-buildx"])
        .run_commands("mkdir /build")
    )
    return image


def main():
    print("Looking up modal.Sandbox app")
    app = modal.App.lookup("docker-test", create_if_missing=True)
    print("Creating sandbox")

    with modal.enable_output():
        sb = modal.Sandbox.create(
            "/usr/bin/dockerd",
            "-D",
            timeout=60 * 60,
            app=app,
            image=create_modal_sandbox_image(),
            experimental_options={"vm_runtime": True},
        )

    print(f"sandbox_id: {sb.object_id}")
    task_id = sb._get_task_id()
    print(f"task_id: {task_id}")
    print(f"To shell into the task, run: modal shell {task_id}")
    # dockerd is the sandbox entrypoint and takes a moment to bind
    # /var/run/docker.sock after the sandbox is created. Poll until the
    # daemon answers so the first `docker build` doesn't run before dockerd is ready.
    print("Waiting for dockerd to be ready")
    wait_p = sb.exec(
        "sh",
        "-c",
        "for i in $(seq 1 120); do "
        "if [ -S /var/run/docker.sock ] && docker info >/dev/null 2>&1; then "
        "echo ready; exit 0; fi; sleep 1; done; "
        "echo 'dockerd not ready after 120s' >&2; exit 1",
    )
    wait_p.wait()
    if wait_p.returncode != 0:
        raise Exception(f"dockerd never became ready: {wait_p.stderr.read()}")

    # A simple Dockerfile that we'll build and run within Modal.
    dockerfile = """
    FROM ubuntu
    RUN apt-get update
    RUN apt-get install -y cowsay curl
    RUN mkdir -p /usr/share/cowsay/cows/
    RUN curl -o /usr/share/cowsay/cows/docker.cow https://raw.githubusercontent.com/docker/whalesay/master/docker.cow
    ENTRYPOINT ["/usr/games/cowsay", "-f", "docker.cow"]
    """
    sb.filesystem.write_text(dockerfile, "/build/Dockerfile")

    print("Building docker image")
    p = sb.exec("docker", "build", "-t", "whalesay", "/build")
    for l in p.stdout:
        print(l, end="")
    p.wait()
    print("--------------------------------")
    if p.returncode != 0:
        print(p.stderr.read())
        raise Exception("Docker build failed")

    # The Sandbox will run a container from the built image and print this:
    #
    #  ________
    # < Hello! >
    #  --------
    #     \
    #      \
    #       \
    #                     ##         .
    #               ## ## ##        ==
    #            ## ## ## ## ##    ===
    #        /"""""""""""""""""\___/ ===
    #       {                       /  ===-
    #        \______ O           __/
    #          \    \         __/
    #           \____\_______/

    print("Running Docker image")
    # Note we can't use -it here because we're not in a TTY.
    p = sb.exec("docker", "run", "--rm", "whalesay", "Hello!")
    print(p.stdout.read())
    p.wait()
    if p.returncode != 0:
        raise Exception(f"Docker run failed: {p.stderr.read()}")
    sb.terminate()

if __name__ == "__main__":
    main()
```

<!-- synmon-sync:docker_in_modal:end -->

</Collapsible>

Additionally, quickly provision a VM Sandbox with a PTY shell via the CLI using:

```
modal shell --experimental-option vm_runtime=1
```

## Running custom init systems

By default, Modal runs and manages the init process (PID 1) inside the VM. Set the `vm_init`
experimental option to an absolute path (typically, `/sbin/init`) to run a conventional init
system such as [`systemd`](https://man7.org/linux/man-pages/man1/systemd.1.html) or
[`openrc`](https://github.com/OpenRC/openrc) as PID 1 instead. Modal's agent runs alongside the
`vm_init`-specified init process.

```python fixture:sb_app
image = modal.Image.from_registry("debian:bookworm-slim").dockerfile_commands(
    "RUN apt-get update",
    # udev is needed so that getty doesn't block for 90s waiting on device activation.
    "RUN apt-get install -y systemd systemd-sysv dbus udev",
    # An empty machine-id tells systemd to generate one at first boot.
    # systemd won't start if the file is missing entirely.
    "RUN rm -f /etc/machine-id && touch /etc/machine-id",
)

sb = modal.Sandbox.create(
    app=sb_app,
    image=image,
    cpu=4,
    memory=2048,
    readiness_probe=modal.Probe.with_exec(
        "systemctl", "is-system-running", interval_ms=250
    ),
    experimental_options={"vm_runtime": True, "vm_init": "/sbin/init"},
)

try:
    sb.wait_until_ready()
finally:
    sb.terminate()
```

## Improvements over gVisor sandboxes

Docker workloads behave more like they do in a non-container environment. In particular:

* Docker state (e.g. `/var/lib/docker`) is included in [Filesystem Snapshots](/docs/guide/sandbox-snapshots#filesystem-snapshots)
* Docker features that previously needed special treatment on gVisor (e.g. inter-container networking) will also work normally

Features that only make sense in a bona fide Linux environment are now available:

* Custom [init systems](https://arxiv.org/pdf/0706.2748) (such as [`systemd`](https://man7.org/linux/man-pages/man1/systemd.1.html)) are supported
* [eBPF](https://ebpf.io/) is supported
* [FUSE](https://www.kernel.org/doc/html/latest/filesystems/fuse.html) mounts are supported
* Resource isolation within the Sandbox via [cgroups](https://man7.org/linux/man-pages/man7/cgroups.7.html) is supported

Finally, for most workloads, the root filesystem will perform better on a VM Sandbox than in a gVisor Sandbox.

## Resource model

Unlike [resource provisioning](/docs/guide/resources) in other runtimes,
memory provisioning is **static** for VM Sandboxes: you get exactly as much
RAM as you request via `memory` argument to `Sandbox.create`. By default, VM
sandboxes get 1GiB of RAM.

However, CPU provisioning is elastic. You can burst above your requested amount.

Costs for both resources are calculated based on the requested amount, used amount,
the duration of Sandbox execution, and [our rates for `cpu` and `memory`](/pricing).

## Limitations

The following limitations are known and we're tracking them:

* **GPUs are not supported.** VM Sandboxes currently only support CPU workloads.
* **The [Sandbox filesystem API](/docs/guide/sandbox-files#filesystem-api-beta) is only available in new SDK versions**. For the Python SDK, it requires version ≥ 1.4.0 and for the JS/TS/Go SDKs, it requires versions ≥ 0.7.6.
* **[`Sandbox.reload_volumes()`](/docs/sdk/py/latest/Sandbox#reload_volumes) is not supported.** VM Sandboxes do not currently support reloading volumes at runtime.
* **Nested KVM virtualization is disabled by default.** Please reach out to us if you'd like access.
* **[VM Memory Snapshots](/docs/guide/vm-memory-snapshots) are only available to a set of enabled customers.**
  Please reach out to us if you'd like early access.
* **Root images ≥ 512 GiB are not supported.** The VM root filesystem is currently limited to 512 GiB. Sandboxes created from container images exceeding this size will fail to start.

If you hit a rough edge that isn't listed here, please reach out via [Slack](/slack) or email us at <support@modal.com>.
