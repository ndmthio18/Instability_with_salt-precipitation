# 本次对话文件归档（部分完成）

归档分支：`figures/prat-interface-editable`。只新增归档和还原工具，不修改默认分支、模型代码或计算结果。

## 优先使用的可编辑文件

| 文件 | 内容 |
|---|---|
| [prat_interface_editable.svg](editable/prat_interface_editable.svg) | 最新土壤–盐壳界面通量与饱和状态一致性诊断图，分层可编辑 |
| [Prat_crust_salt_balance_editable.svg](editable/Prat_crust_salt_balance_editable.svg) | 盐质量收支示意图的原始可编辑版本 |
| [prat_editable_direct.svg](editable/prat_editable_direct.svg) | 在对话中直接提供源码的精简排版版本 |
| [prat_salt_balance.svg](editable/prat_salt_balance.svg) | 原始收支图的备用文件名，内容未改动 |
| [prat_interface_source.txt](editable/prat_interface_source.txt) | 最新界面图的完整 SVG 源码 |
| [prat_editable_direct.txt](editable/prat_editable_direct.txt) | 精简排版版本的完整 SVG 源码 |
| [prat_salt_balance.zip](original_archives/prat_salt_balance.zip) | 与对话中小型 ZIP 的 SHA-256 完全一致 |
| [figure9_physical_columns.png](images/figure9_physical_columns.png) | 对话上传的计算图；复用仓库已有同一 Git blob |

[下载已核验文件的新 ZIP 包](verified_download_bundle.zip) · [核验清单](VERIFIED_FILES.json) · [全部文件上传状态](UPLOAD_STATUS.md)

**这不是全部对话附件的完整备份。** 未上传的原始 PNG、预览和大 ZIP 明确列在上传状态表中，不能把再渲染的图片当作这些原始附件。

## 如何打开

在 GitHub 文件页面使用 Download raw file；SVG 在浏览器中打开后也可另存为。SVG 内的文字和矢量线条保持原样，没有将文字转成轮廓。可使用矢量编辑器重新排版。历史图件包含讨论过程中的近似与条件，不应将归档本身视为科学验证。

## 完整性与还原

六份 SVG/TXT 原始字节通过无损压缩 payload 传输。脚本在写文件前验证压缩包及每个文件的 SHA-256，只写当前归档目录。可在本仓库运行：

```bash
python3 conversation_files/2026-11-10/restore_editables.py
```

如果 `editable/` 尚未出现，说明分支内的还原作业尚未成功；源文件字节已保存在 `_payload/` 中，不需要访问 ChatGPT 的临时目录。
