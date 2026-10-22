# Service Users

<Callout variant="gated-feature">
Service users are available on the <a href="/pricing">Team and Enterprise plans</a>. Visit <a href="/settings/plans">Workspace settings</a> to upgrade.
</Callout>

Service users are programmatic accounts that allow automated systems to interact with Modal. They're ideal for CI/CD pipelines, automated deployments, and other workflows that need to authenticate.

## Create a Service User

Service users are only available for shared workspaces. You will need workspace owner or manager privileges to create service users.

To create a service user:

1. Go to your workspace [tokens settings page](/settings/tokens/service-users)
2. Click **New Service User**
3. Enter a name for your service user (must be lowercase alphanumeric, can contain hyphens or underscores)
4. Click **Create**
5. Copy the `MODAL_TOKEN_ID` and `MODAL_TOKEN_SECRET`. **This is the only time you can view the token secret** for security reasons.
6. Click **Configure environments** to grant the service user access. A new service user has **No Access** to every Environment, so it can't be used until you assign it the **Viewer** or **Contributor** role on at least one Environment.

## Use Service User Tokens

Set the service user credentials as environment variables in your automated environment:

```bash
export MODAL_TOKEN_ID=your-token-id
export MODAL_TOKEN_SECRET=your-token-secret
```

Once configured, you can use Modal's CLI and Python SDK as usual:

```bash
modal deploy your_app.py
```

## Delete a Service User

To remove a service user:

1. Go to the [tokens settings page](/settings/tokens/service-users)
2. Find the service user in the table
3. Click **Delete** when you hover over the row

## Permissions

Service users default to **No Access** on every Environment — they cannot read or write to an Environment until you grant them a Role. This differs from Workspace Members, who default to **Contributor**.

Grant a service user the **Viewer** or **Contributor** role on the specific Environments it needs, either:

* from the service user: on the [tokens settings page](/settings/tokens/service-users), open the **⋯** menu on the service user's row and select **Manage environments**, or
* from an Environment: on the [Environment settings page](/settings/workspace-management/environments), click **Manage** on the Environment, then open the **Access Restrictions** tab.

See [Role-Based Access Control](/docs/guide/rbac) for what each Environment Role allows.

## Securing Service Users

Because service user tokens are long-lived and used in automated environments, it's important to limit their access to only what's necessary:

* **Store tokens securely.** Use a secrets manager or your CI/CD platform's built-in secrets storage rather than hardcoding tokens in source code or configuration files.
* **Grant least-privilege access.** Because service users start with **No Access**, assign the **Viewer** or **Contributor** role only on the Environments they actually need, keeping production isolated from development and staging.
