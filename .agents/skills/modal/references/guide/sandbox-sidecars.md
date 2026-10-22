# Sandbox Sidecars

<Callout variant="alpha">

There are currently several [known limitations](#limitations).

</Callout>

## Introduction

Sandbox Sidecars let you run additional containers alongside your main
Sandbox container, on the same host. A sandbox and its sidecars are connected
via an internal bridge network, allowing low latency communication between
containers over TCP/UDP, making them ideal for:

* Separating an agent harness from its execution environment, by running the
  agent in one container and its tool calls in another
* Credentials injection, by running a proxy in a separate, trusted container
  from the primary application, and letting that proxy inject credentials or
  other secrets before passing on network calls to external services.
  See the [secrets injection example](/docs/examples/sidecar_secrets_injection)
  for a working demonstration
* Splitting out complex multi-service applications over separate containers,
  such as databases, caches or worker processes, similar to Docker Compose.

We're still discovering all the ways that Sandbox Sidecars can be used - if you
come up with another use case, please let us know!

Sidecars are managed through the sidecars interface on a Sandbox
(`_experimental_sidecars` in Python, `experimentalSidecars` in JS/Go),
which provides methods to create, list, get, and terminate Sidecar containers.

Each Sidecar container:

* Runs its own image independently from the main Sandbox container.
* Runs in a separate, sandboxed process isolated from the main Sandbox container and other Sidecar containers.
* Can communicate over an internal bridge network with the main Sandbox container and other Sidecar containers.
* Can be created, terminated, and replaced dynamically during the Sandbox's lifetime.
* Supports executing commands just like the main Sandbox container.

## Usage

### Creating a Sidecar container

The main Sandbox container is resolvable as `main`, and each Sidecar container
is resolvable by the `name` you give it at creation time.

<CodeTabs>
{#snippet python()}

```python notest
import modal

app = modal.App.lookup("sidecar-example", create_if_missing=True)
image = modal.Image.debian_slim().build(app)

sb = modal.Sandbox.create("sleep", "600", app=app, image=image, timeout=300)

sidecar = sb._experimental_sidecars.create(
    "python",
    "-m",
    "http.server",
    "8080",
    name="web",
    image=image,
)

# Give the server a moment to start, then call it from the main sandbox.
p = sb.exec(
    "python",
    "-c",
    "import time, urllib.request; time.sleep(1); print(urllib.request.urlopen('http://web:8080').status)",
)
p.wait()
print(p.stdout.read())  # "200"

sb.terminate()
```

{/snippet}

{#snippet javascript()}

```javascript notest
import { ModalClient } from "modal";

const modal = new ModalClient();
const app = await modal.apps.fromName("sidecar-example", {
  createIfMissing: true,
});
const image = await modal.images.fromRegistry("python:3.13-slim").build(app);

const sb = await modal.sandboxes.create(app, image, {
  command: ["sleep", "600"],
  timeoutMs: 300 * 1000,
});

const sidecar = await sb.experimentalSidecars.create("web", image, {
  command: ["python", "-m", "http.server", "8080"],
});

// Give the server a moment to start, then call it from the main sandbox.
const p = await sb.exec([
  "python",
  "-c",
  "import time, urllib.request; time.sleep(1); print(urllib.request.urlopen('http://web:8080').status)",
]);
await p.wait();
console.log(await p.stdout.readText()); // "200"

await sb.terminate();
```

{/snippet}

{#snippet go()}

```go
package main

import (
	"context"
	"fmt"
	"io"
	"time"

	modal "github.com/modal-labs/modal-client/go"
)

func main() {
	ctx := context.Background()
	mc, _ := modal.NewClient()

	app, _ := mc.Apps.FromName(ctx, "sidecar-example", &modal.AppFromNameParams{
		CreateIfMissing: true,
	})
	image, _ := mc.Images.FromRegistry("python:3.13-slim", nil).Build(ctx, app, nil)

	sb, _ := mc.Sandboxes.Create(ctx, app, image, &modal.SandboxCreateParams{
		Command: []string{"sleep", "600"},
		Timeout: 5 * time.Minute,
	})
	defer sb.Terminate(ctx, nil)

	sidecar, _ := sb.ExperimentalSidecars.Create(ctx, "web", image, &modal.SidecarCreateParams{
		Command: []string{"python", "-m", "http.server", "8080"},
	})
	_ = sidecar

	// Give the server a moment to start, then call it from the main sandbox.
	p, _ := sb.Exec(ctx, []string{
		"python", "-c",
		"import time, urllib.request; time.sleep(1); print(urllib.request.urlopen('http://web:8080').status)",
	}, nil)
	stdout, _ := io.ReadAll(p.Stdout)
	fmt.Println(string(stdout)) // "200"
}
```

{/snippet} </CodeTabs>

Names are resolved using `/etc/hosts` which gets updated when a sidecar is created or terminated.

### Listing and retrieving sidecars

You can list all running Sidecar containers or retrieve a specific one by name:

<CodeTabs>
{#snippet python()}

```python notest
containers = sb._experimental_sidecars.list()
for container in containers:
    print(f"{container.name}: {container.object_id}")

sidecar = sb._experimental_sidecars.get(name="web")
```

{/snippet}

{#snippet javascript()}

```javascript notest
const containers = await sb.experimentalSidecars.list();
for (const container of containers) {
  console.log(`${container.containerName}: ${container.containerId}`);
}

const sidecar = await sb.experimentalSidecars.get("web");
```

{/snippet}

{#snippet go()}

```go notest
containers, _ := sb.ExperimentalSidecars.List(ctx, nil)
for _, container := range containers {
	fmt.Printf("%s: %s\n", container.ContainerName, container.ContainerID)
}

sidecar, _ := sb.ExperimentalSidecars.Get(ctx, "web", nil)
_ = sidecar
```

{/snippet} </CodeTabs>

### Routing HTTPS traffic through a Sidecar

Sidecars can be used to inspect the outgoing HTTPS traffic from the main Sandbox
container in a different context, for example to perform more advanced request
filtering, inspecting requests for logging, or injecting secrets that the main
Sandbox container should not have access to.

Normally, applications need to support explicit proxy configuration, such as
respecting the `HTTPS_PROXY` environment variable, to route traffic through a Sidecar.
To include HTTPS traffic (TCP on port 443) from **proxy-unaware** applications, you can
set the experimental `proxy_traffic_via_sidecar` option that routes **all** outbound
HTTPS traffic from the main Sandbox container through a Sidecar.

<CodeTabs>
{#snippet python()}

```python notest
sb = modal.Sandbox.create(
    "sleep",
    "600",
    app=app,
    image=image,
    experimental_options={"proxy_traffic_via_sidecar": "my-proxy-sidecar"},
)

# Until this Sidecar is running, HTTPS from the main container is refused.
sb._experimental_sidecars.create("python", "/proxy.py", name="my-proxy-sidecar", image=proxy_image)
```

{/snippet}

{#snippet javascript()}

```javascript notest
const sb = await modal.sandboxes.create(app, image, {
  command: ["sleep", "600"],
  experimentalOptions: { proxy_traffic_via_sidecar: "my-proxy-sidecar" },
});

// Until this Sidecar is running, HTTPS from the main container is refused.
await sb.experimentalSidecars.create("my-proxy-sidecar", proxyImage, {
  command: ["python", "/proxy.py"],
});
```

{/snippet}

{#snippet go()}

```go notest
sb, _ := mc.Sandboxes.Create(ctx, app, image, &modal.SandboxCreateParams{
	Command:             []string{"sleep", "600"},
	ExperimentalOptions: map[string]any{"proxy_traffic_via_sidecar": "my-proxy-sidecar"},
})
// Until this Sidecar is running, HTTPS from the main container is refused.
sb.ExperimentalSidecars.Create(ctx, "my-proxy-sidecar", proxyImage, &modal.SidecarCreateParams{
	Command: []string{"python", "/proxy.py"},
})
```

{/snippet} </CodeTabs>

The Sidecar receives a raw TLS stream and must read the destination hostname
from the `ClientHello`'s SNI. The original destination IP is not forwarded, so
mechanisms such as `SO_ORIGINAL_DST` do not work. To read or rewrite HTTP
requests, terminate TLS in the proxy using a certificate authority that the
Sandbox trusts. See the [Sidecar traffic routing
example](/docs/examples/sidecar_traffic_routing) for a complete mitmproxy-based
request filter.

Only TCP traffic to port 443 is relayed. The Sandbox's own egress controls are
set aside for it: an `outbound_cidr_allowlist` on the Sandbox still governs every
other port, but relayed traffic passes regardless of what it lists. Non-relayed
traffic is still subject to the egress controls of the Sandbox.

Relayed traffic is instead governed by the egress controls of the Sidecar it is
relayed into. A Sidecar's outbound network policy is independent of the main
container's and defaults to open, so unless you pass `outbound_cidr_allowlist`
or `outbound_domain_allowlist` to the Sidecar itself, relayed traffic reaches
any destination the Sidecar chooses to connect to.

The option cannot be combined with setting `block_network`,
`outbound_domain_allowlist` or `proxy` on the Sandbox.

### Filesystem snapshots

You can snapshot a running Sidecar's filesystem into a reusable Image. The
resulting Image can be used anywhere an existing Image is accepted; the example
below uses it to start another Sidecar. The snapshot is scoped to that Sidecar;
it does not include the main Sandbox filesystem or other Sidecars.

<CodeTabs>
{#snippet python()}

```python notest
sidecar.filesystem.write_text("ready", "/tmp/state")
snapshot = sidecar.snapshot_filesystem()

restored = sb._experimental_sidecars.create(
    "sleep", "600", name="restored", image=snapshot
)
assert restored.filesystem.read_text("/tmp/state") == "ready"
```

{/snippet}

{#snippet javascript()}

```javascript notest
await sidecar.filesystem.writeText("ready", "/tmp/state");
const snapshot = await sidecar.snapshotFilesystem();

const restored = await sb.experimentalSidecars.create("restored", snapshot, {
  command: ["sleep", "600"],
});
console.assert((await restored.filesystem.readText("/tmp/state")) === "ready");
```

{/snippet}

{#snippet go()}

```go notest
_ = sidecar.Filesystem.WriteText(ctx, "ready", "/tmp/state", nil)
snapshot, _ := sidecar.SnapshotFilesystem(ctx, nil)

restored, _ := sb.ExperimentalSidecars.Create(ctx, "restored", snapshot, &modal.SidecarCreateParams{
	Command: []string{"sleep", "600"},
})
state, _ := restored.Filesystem.ReadText(ctx, "/tmp/state", nil)
fmt.Println(state) // "ready"
```

{/snippet} </CodeTabs>

## Resource configuration

The main Sandbox container and the Sidecar containers share the resource allocation (CPU and memory) of the Sandbox,
and resources are configured only on the Sandbox. When planning your
resource allocation, make sure the Sandbox is configured with enough CPU
and memory for all containers combined.
Bursting is still possible, see the [guide to Sandbox resources and
pricing](/docs/guide/sandbox-resources) for more details.

For example, if you want to run a Sandbox with two Sidecars, and you expect the main
container to use 1 CPU core and 512 MiB of memory, Sidecar A to use 0.5 CPU and 256 MiB,
and Sidecar B to use 0.5 CPU and 256 MiB, you should set the Sandbox's resources to at
least 2 CPUs and 1024 MiB to accommodate all three containers.

The maximum number of Sidecars you can create is also determined by the main Sandbox's
resource reservation. Each container (including the main one) requires a minimum of
32 mCPU and 32 MiB of memory, so the limit is:

```
max containers = min(cpu_in_milli / 32, memory_in_mib / 32)
```

There is also a hard limit of **250** concurrent sidecar containers per sandbox,
regardless of the resource reservation.

## Limitations

The main sandbox supports the same features as a regular sandbox, but some features are not yet supported
for sidecars:

* **Pre-built images only**: Sidecar images must be pre-built using `image.build()`, referenced
  by ID via `Image.from_id()` or name via `Image.from_name()`, or created from filesystem/directory snapshots. Lazy image
  building is not supported for sidecars. See also [Separating Image builds from Sandbox creation](/docs/guide/sandboxes#separating-image-builds-from-sandbox-creation).
* **No GPU support**: Sidecar containers cannot access GPUs, even when the Sandbox is configured with one.
* **No Cloud Bucket Mount support**: Sidecar containers do not currently support attaching [Cloud Bucket Mounts](/docs/guide/cloud-bucket-mounts).
* **No memory snapshot support**: A Sidecar's filesystem can be snapshotted
  independently, but Sidecar memory state is not captured in
  [Sandbox snapshots](/docs/guide/sandbox-snapshots).
* **VM incompatibility**: Sidecars are not compatible with VM Sandboxes.
* **Changes to /etc/hosts are not preserved**: `/etc/hosts` is rewritten on sidecar create/terminate and user changes are not preserved.
* **Maximum of 250 concurrent sidecars**: A sandbox can have at most 250 sidecar containers running at the same time.
* **No [Proxy](/docs/guide/proxy-ips) support**: Traffic from a Sidecar does not exit through a Proxy. Because relayed traffic leaves from the Sidecar, a Sandbox cannot currently combine a Proxy with `proxy_traffic_via_sidecar`.
