# Proxies

<Callout variant="beta" />

You can securely connect with resources in your private network
using a Modal Proxy. Proxies are a secure tunnel between
Apps and exit nodes with static IPs. You can allow-list those static IPs
in your network firewall, making sure that only traffic originating from these
IP addresses is allowed into your network.

Proxies are unique and not shared between workspaces. All traffic
between your Apps and the Proxy server is encrypted using
[WireGuard](https://www.wireguard.com/).

Modal Proxies are built on top of [vprox](https://github.com/modal-labs/vprox),
a Modal open-source project used to create highly available proxy servers
using WireGuard.

## Creating a Proxy

<Callout variant="gated-feature">

Proxies are available on the <a href="/pricing">Team and Enterprise plans</a>. Visit <a href="/settings/plans">workspace settings</a> to upgrade.

</Callout>

You can create Proxies in your workspace [Settings](/settings) page.
Team Plan users can create one Proxy and Enterprise users three Proxies. Each Proxy
can have a maximum of five static IP addresses.

Please reach out to <support@modal.com> if you need greater limits.

You can select which region the Proxy should live in when creating it.
This allows you to control where in the world the traffic will be routed through.
Placing it closer to the resources you want to connect with can reduce overall latency.

## Using a Proxy

After a Proxy is online, add it to a Modal Function with the argument
`proxy=Proxy.from_name("<your-proxy>")`. For example:

```python
import modal
import subprocess

app = modal.App(image=modal.Image.debian_slim().apt_install("curl"))

@app.function(proxy=modal.Proxy.from_name("<your-proxy>"))
def my_ip():
    subprocess.run(["curl", "-s", "ifconfig.me"])

@app.local_entrypoint()
def main():
    my_ip.remote()
```

All network traffic from your Function will now use the Proxy as a tunnel.

The program above will always print the same IP address independent
of where it runs in Modal's infrastructure. If that same program
were to run without a Proxy, it would print a different IP
address depending on where it runs.

## Proxy performance

All traffic that goes through a Proxy is encrypted by WireGuard. This adds
latency to your Function's networking. If you are experiencing networking issues
with Proxies related to performance, first add more IP addresses to your
Proxy (see [Adding more IP addresses to a Proxy](#adding-more-ip-addresses-to-a-proxy)).

Placing your Proxy geographically closer to a target may also improve overall performance.
See [Creating a Proxy](#creating-a-proxy) for selecting the Proxy region when creating it.

## Adding more IP addresses to a Proxy

Proxies support up to five static IP addresses. Adding IP addresses improves
throughput linearly.

You can add an IP address to your workspace in [Settings](/settings) > Proxies.
Select the desired Proxy and add a new IP.

If a Proxy has multiple IPs, Modal will randomly pick one when running your Function.

## Proxies and Sandboxes

Proxies can also be used with [Sandboxes](/docs/guide/sandboxes). For example:

```python notest
import modal

app = modal.App.lookup("sandbox-proxy", create_if_missing=True)
sb = modal.Sandbox.create(
    app=app,
    image=modal.Image.debian_slim().apt_install("curl"),
    proxy=modal.Proxy.from_name("<your-proxy>"))

process = sb.exec("curl", "-s", "https://ifconfig.me")
stdout = process.stdout.read()
print(stdout)

sb.terminate()
```

Similarly to our Function implementation, this Sandbox program will
always print the same IP address.

A Proxy covers the Sandbox's main container. Traffic from a
[Sidecar](/docs/guide/sandbox-sidecars) does not currently exit through it, so a
Proxy cannot be combined with the experimental `proxy_traffic_via_sidecar`
option, which would relay the main container's HTTPS out through a Sidecar.
