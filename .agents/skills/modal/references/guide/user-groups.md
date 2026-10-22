# User groups

<Callout variant="gated-feature">
User groups require <a href="/docs/guide/rbac">Role-Based Access Control
(RBAC)</a>, which is available on the <a href="/pricing">Team and Enterprise
plans</a>. Creating and managing groups requires the Workspace Owner or
Workspace Manager Role. Visit <a href="/settings/plans">Workspace settings</a>
to upgrade.
</Callout>

<Callout variant="beta" />

User groups let you manage access for multiple Workspace Members at once. Add
members to a group, then assign the group an Environment Role in each
[Restricted Environment](/docs/guide/rbac#setting-up-restricted-environments)
it needs to access.

## Group types

Modal supports two types of user groups:

| Group type | Created in | Name and membership managed in |
| ---------- | ---------- | ------------------------------ |
| Custom     | Modal      | Modal                          |
| SCIM       | Your IdP   | Your IdP                       |

You can assign Environment Roles to either type of group from Modal. SCIM only
controls the group's name and membership; Modal remains the source of truth for
its Environment access.

## Create a custom group

To create and populate a custom group:

1. Open the **Groups** tab on the [Workspace Management settings
   page](/settings/workspace-management/groups).
2. Click **Create Group**.
3. Enter a group name and select the Workspace Members to add.
4. Click **Create**.

A group can be created without any members or Environment access. Assigning an
Environment Role later applies that Role to every current and future member of
the group.

## Update a custom group

On the **Groups** settings page, expand a custom group and click **Edit Group**
to rename it or change its membership. You can also remove a member directly
from the group's member list.

Changes to group membership update the members' access to every Environment
assigned to that group.

## Manage groups with SCIM

When [SCIM](/docs/guide/scim) group provisioning is enabled, groups pushed by
your identity provider (IdP) appear automatically on the **Groups** settings
page with a **SCIM** label.

Manage the name and membership of a SCIM group in the IdP. These fields cannot
be edited in Modal, so the IdP remains their source of truth. You can still
assign, change, and remove the group's Environment Roles in Modal.

Removing a user from a group in the IdP removes the access they received from
that group after the change syncs to Modal. Deactivating or removing the user
through SCIM also removes their Workspace access according to your SCIM
configuration.

## Assign an Environment Role to a group

User groups work with [Role-Based Access Control
(RBAC)](/docs/guide/rbac) at the Environment level. Groups do not receive a
Workspace Role and do not change who can manage Workspace settings or billing.

To assign a group to an Environment:

1. Open [Environment settings](/settings/workspace-management/environments).
2. Click **Manage** for a Restricted Environment, then open **Access
   Restrictions**.
3. Select the **Groups** tab and click **Add Group**.
4. Select a custom or SCIM group and choose an Environment Role:
   * **Contributor** — can view and modify resources in the Environment.
   * **Viewer** — can view resources but cannot modify them.
   * **No Access** — cannot discover or access the Environment.
5. Click **Confirm**.

Group Roles can only be assigned in Restricted Environments. From the same
**Groups** tab, you can change a group's Role or remove the group from the
Environment. The **Groups** settings page also lists every Environment and Role
assigned to each group.

## How group Roles combine with other RBAC settings

A Workspace Member can receive an Environment Role directly, through one or
more groups, or from the Restricted Environment's default member Role. Modal
resolves these sources as follows:

* Any direct or group Role takes precedence over the default member Role.
* If a member has multiple direct or group Roles, the most permissive Role
  applies: **Contributor**, then **Viewer**, then **No Access**.
* Workspace Owners and Managers always have **Contributor** access, regardless
  of group assignments.

For example, a member assigned **Viewer** directly and **Contributor** through a
group has Contributor access. A **No Access** group assignment does not revoke a
more permissive direct or group assignment, but it does override an Environment
default of Viewer or Contributor when no other explicit assignment applies.

To keep access predictable, assign each team to the smallest set of groups it
needs and review both direct and group Roles when troubleshooting a member's
effective access.
