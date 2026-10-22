# SCIM Integration

<Callout variant="gated-feature">
SCIM support is available on the <a href="/pricing">Enterprise plan</a>. Contact <a href="mailto:sales@modal.com">sales@modal.com</a> for more information.
</Callout>

<Callout variant="beta" />

[SCIM (System for Cross-domain Identity Management)](https://datatracker.ietf.org/doc/html/rfc7643) is a protocol that Identity Providers (IdPs) can use to automate user management in connected apps.

Modal supports SCIM for automatic provisioning and deprovisioning of users.

## Connecting an IdP

### Step 1: Generate a SCIM token

1. Sign in to https://modal.com and visit your [Workspace Management](/settings/workspace-management/identity-and-provisioning) page's "Identity and Provisioning" tab. If SCIM is enabled for your Workspace, there will be a "SCIM Tokens" section on the page below the SSO configuration settings. If you do not see a section for SCIM tokens, contact Modal support about enabling SCIM support for your Workspace.
2. Click on "New SCIM Token" then "Create Token".
3. A new token will be generated and displayed to you. Copy the value from the "Token Secret" box and store it somewhere secure. You can also copy the exact url that the IdP will require to integrate with your Modal Workspace. Once you click "Done", you will not be able to view the token secret again and will have to generate a new one if you can't otherwise access it.

### Step 2: IdP Configuration

#### Okta

1. Create a private SCIM integration.

   The SCIM integration must be separate from the Modal catalog app used for
   SAML SSO. Your existing Modal app can continue to handle SSO; do not add SSO
   to the SCIM integration.

   In the Okta Admin Console:

   1. Go to "Applications > Applications".
   2. Click "Create a new app integration".
   3. Select "Okta Integration Wizard".
   4. Choose "Provisioning" as the capability.
   5. Choose "SCIM 2.0" as the provisioning method.

2. Configure the integration with your Modal SCIM credentials.

   In the wizard's provisioning settings, enter:

   | Okta setting                      | Modal value                                           |
   | --------------------------------- | ----------------------------------------------------- |
   | SCIM connector base URL           | `https://modal.com/api/<your-workspace>/scim/v2`      |
   | Unique identifier field for users | `userName`                                            |
   | Authentication mode               | HTTP Header                                           |
   | Authorization                     | The full SCIM token generated in step 1               |
   | Supported provisioning actions    | Push New Users, Push Profile Updates, and Push Groups |

   Test the API credentials, then review and deploy the integration. When
   prompted, add an app instance from your organization's "Private apps"
   catalog.

   In the app instance's "Provisioning > To App" settings, enable "Create
   Users", "Update User Attributes", and "Deactivate Users". Assign the
   people and groups that Okta should provision to Modal.

   For more information, see Okta's
   [Okta Integration Wizard documentation](https://help.okta.com/en-us/Content/Topics/Apps/oiw/create-app-integration.htm).

#### Microsoft Entra ID

1. Create or select an Enterprise Application.

   In the [Microsoft Entra admin center](https://entra.microsoft.com/):

   * If you already have a Modal Enterprise Application that you use for SSO,
     you can reuse it for SCIM provisioning. Go to "Entra ID > Enterprise apps"
     and select the application.
   * Otherwise, go to "Entra ID > Enterprise apps" and create an application:
     1. Select "New application > Create your own application".
     2. Enter a name such as `Modal SCIM`.
     3. Select "Integrate any other application you don't find in the gallery
        (Non-gallery)" and create the application.

2. Configure the application with your Modal SCIM credentials.

   Open the Enterprise Application, select "Provisioning > New configuration",
   and enter:

   | Entra setting | Modal value                                      |
   | ------------- | ------------------------------------------------ |
   | Tenant URL    | `https://modal.com/api/<your-workspace>/scim/v2` |
   | Secret Token  | The full SCIM token generated in step 1          |

   Select "Test Connection", then create the provisioning configuration.
   Review the user attribute mappings and ensure that an email address is mapped
   to `userName`. Assign the users and groups that Entra should provision, then
   select "Start provisioning".

   For more information, see Microsoft's
   [SCIM provisioning documentation](https://learn.microsoft.com/en-us/entra/identity/app-provisioning/use-scim-to-provision-users-and-groups#integrate-your-scim-endpoint-with-the-microsoft-entra-provisioning-service).

#### Other IdPs

1. Create a custom SCIM integration.

   Create a custom or non-gallery application that supports outbound SCIM 2.0
   provisioning. Name it something recognizable, such as `Modal SCIM`. An
   existing SSO integration can continue to handle authentication; whether the
   SCIM integration must be a separate app depends on your IdP.

2. Configure the integration with your Modal SCIM credentials.

   Find your IdP's equivalent settings and enter:

   | Setting                | Modal value                                      |
   | ---------------------- | ------------------------------------------------ |
   | SCIM version           | 2.0                                              |
   | SCIM base URL          | `https://modal.com/api/<your-workspace>/scim/v2` |
   | Authorization method   | Bearer token                                     |
   | Token                  | The full SCIM token generated in step 1          |
   | Unique user identifier | Email address in `userName`                      |

   Enable creating, updating, and deactivating users. You can also enable group
   provisioning. Test the connection, assign the users and groups that your IdP
   should provision, and start provisioning.

Modal supports the following SCIM capabilities:

| Capability                         | Supported | Notes                         |
| ---------------------------------- | --------- | ----------------------------- |
| SCIM 2.0                           | Yes       |                               |
| Pagination                         | Yes       |                               |
| Create, update, and remove users   | Yes       |                               |
| Create, update, and delete groups  | Yes       |                               |
| Update group membership with PATCH | Yes       |                               |
| Generate temporary passwords       | No        | Modal authentication uses SSO |

The IdP may also ask you to specify which user attributes are supported.

| SCIM user attribute | Modal support | Notes                                           |
| ------------------- | ------------- | ----------------------------------------------- |
| externalId          | Yes           |                                                 |
| userName            | Yes           | Required; must contain the user's email address |
| displayName         | Yes           |                                                 |
| name.familyName     | Yes           |                                                 |
| name.givenName      | Yes           |                                                 |
| emails              | Read only     | The primary email is derived from `userName`    |
| active              | Yes           |                                                 |
| addresses           | No            |                                                 |
| profileUrl          | No            |                                                 |

## Managing Tokens

Token management is restricted to only workspace owners and managers.

Up to two SCIM tokens may be active at any time. It may be useful to generate a second token to facilitate seamless token rotation - a workspace admin can generate a new token, use it to replace the old one in the connected IdP, and finally revoke the old token to ensure that no updates are dropped during the rotation process. Except during the process of rotation we recommend having only one SCIM token active at a time as a security best practice.

## Troubleshooting

If your IdP indicates that it is unable to authenticate with Modal, first double check that the token was copied correctly - the full token will have the form `si-XXXXXXXXXXXXXXXXXXXXXX:ss-XXXXXXXXXXXXXXXXXXXXXX`.

If you experience any issues with or have any questions about SCIM integration, please reach out via [Slack](/slack) or email us at <support@modal.com>.
