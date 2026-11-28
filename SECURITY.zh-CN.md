# 安全策略

[English](SECURITY.md) | [简体中文](SECURITY.zh-CN.md)

## 支持的版本

| 版本 | 安全更新 |
|---|---|
| `3.0.0-beta.1` | 支持 |

## 报告漏洞

请通过 **security@apich.org** 私下报告疑似漏洞，并注明受影响的 RSSN 版本、组件和影响范围；在可以安全分享的前提下请附最小复现样例。

## 嵌入边界

- `plugin-cdylib` 会在宿主进程中加载本地代码，请仅用于可信的共享库。
- C ABI 要求调用者提供有效的指针、长度和句柄生命周期；公开头文件记录了每个返回值的所有权归属。
- 请在有限的搜索预算和受限的宿主资源下解析并执行不可信输入。如发现异常资源消耗或内存安全缺陷，请通过上述私密渠道报告。

## 组织贡献者标准

以下要求适用于组织成员；外部贡献者可将其作为配置自身开发环境的参考。

- **工作站**：推荐使用 Fedora KDE Plasma Desktop，并在受组织管理的工作站上启用 LUKS 全盘加密。
- **终端安全**：组织主机启用 Wazuh 代理、强制模式 SELinux、firewalld、Lynis、ClamAV、rkhunter、chkrootkit、AIDE 与 unhide，并配合 Wazuh、Greenbone 及经授权的网络资产盘点服务。
- **身份认证**：组织账号与 SSH 访问须使用基于 YubiKey 的身份认证。
