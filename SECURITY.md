# Security Policy

## Supported versions · 支持的版本

| Version · 版本 | Security updates · 安全更新 |
| --- | --- |
| `3.0.0-beta.1` | Supported · 支持 |

## Reporting a vulnerability · 报告漏洞

Report suspected vulnerabilities privately to **security@apich.org**. Describe the affected RSSN version, component and impact, and include a minimal supporting example when it is safe to share.

请通过 **security@apich.org** 私下报告疑似漏洞，并注明 RSSN 版本、受影响的组件和问题影响。可以安全分享时请附最小样例。

## Embedding boundaries · 嵌入边界

- `plugin-cdylib` loads native code into the host process. Enable it for trusted shared libraries.
- The C ABI requires valid pointers, lengths and handle lifetimes. Its public headers document ownership of each return value.
- Parse and execute untrusted input with a finite search budget and bounded host resources. Report unexpected resource use and memory-safety defects through the private channel.

- `plugin-cdylib` 在宿主进程加载本地代码。仅用于可信共享库。
- C ABI 要求调用者提供有效指针、长度和句柄生存期。公开头文件说明每个返回值的所有权。
- 在有限搜索预算和受限宿主资源中解析、执行不可信输入。意外资源消耗和内存安全问题可通过上述私下渠道报告。

## Organization contributor standards · 组织贡献者标准

These requirements apply to organization members. External contributors can use them as guidance for their own development environments.

- **Workstations.** Fedora KDE Plasma Desktop is the recommended environment. Enable LUKS full-disk encryption when installing a managed workstation.
- **Endpoint security.** The organization uses Wazuh agents, enforcing SELinux, firewalld and Lynis, ClamAV, rkhunter and chkrootkit, plus AIDE and unhide for integrity checks. Managed security services include Wazuh, Greenbone and authorized network inventory. Apply host-specific settings and scan schedules through the organization’s operations process.
- **Identity.** Organization accounts and SSH access require YubiKey-backed authentication.

这些要求适用于组织成员。外部贡献者可以参考这些要求配置自己的开发环境。

- **工作站。** 推荐使用 Fedora KDE Plasma Desktop。安装受组织管理的工作站时启用 LUKS 全盘加密。
- **终端安全。** 组织使用 Wazuh 代理、强制模式的 SELinux、firewalld、Lynis、ClamAV、rkhunter 和 chkrootkit，并使用 AIDE 和 unhide 检查系统完整性。组织管理的安全服务包括 Wazuh、Greenbone 和经授权的网络资产盘点。通过组织运维流程配置各主机的设置和扫描计划。
- **身份认证。** 组织账号和 SSH 访问须使用 YubiKey 支持的身份认证。
