# Security Policy

[English](SECURITY.md) | [简体中文](SECURITY.zh-CN.md)

## Supported Versions

| Version | Security updates |
|---|---|
| `3.0.0-beta.1` | Supported |

## Reporting a Vulnerability

Report suspected vulnerabilities privately to **security@apich.org**. Describe the affected RSSN version, component, and impact, and include a minimal supporting example when it is safe to share.

## Embedding Boundaries

- `plugin-cdylib` loads native code into the host process. Enable it only for trusted shared libraries.
- The C ABI requires valid pointers, lengths, and handle lifetimes. Public headers document ownership of each return value.
- Parse and execute untrusted input within a finite search budget and bounded host resources. Report unexpected resource consumption or memory-safety defects through the private channel above.

## Organization Contributor Standards

These requirements apply to organization members; external contributors may use them as guidance for their own environments.

- **Workstations**: Fedora KDE Plasma Desktop is recommended, with LUKS full-disk encryption enabled on managed workstations.
- **Endpoint security**: Managed hosts use Wazuh agents, enforcing SELinux, firewalld, Lynis, ClamAV, rkhunter, chkrootkit, AIDE, and unhide, alongside Wazuh, Greenbone, and authorized network inventory services.
- **Identity**: Organization accounts and SSH access require YubiKey-backed authentication.
