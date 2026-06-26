# Mini Randomized & Derandomization（迷你随机化与去随机化）

**从零开始、零依赖的 C 语言实现**，涵盖大学级随机化计算与去随机化理论。每个模块对应 MIT（及其他顶尖大学）的课程，将随机性提取、伪随机生成和硬度与随机性定理转化为可运行的 C 代码。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|-----------|--------|-------------|
| [mini-bpp-rp-zpp-classes](mini-bpp-rp-zpp-classes/) | 随机化复杂度类（BPP/RP/ZPP/co-RP）、概率放大、Chernoff 界、Schwartz-Zippel 多项式恒等测试、Freivalds 验证 | MIT 6.045, MIT 6.841, Stanford CS254 |
| [mini-derandomizing-space](mini-derandomizing-space/) | 空间界限去随机化、Savitch 定理、Immerman-Szelepcsényi 定理（NL=coNL）、Reingold USTCON、Nisan 空间 PRG、分支程序 | MIT 6.841/18.404, Stanford CS254, Berkeley CS278 |
| [mini-expander-graphs-construction](mini-expander-graphs-construction/) | 谱扩展图、Cheeger 不等式、zigzag 积、Ramanujan 图（LPS）、Margulis/Gabber-Galil 构造、扩展混合引理 | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-hardness-vs-randomness](mini-hardness-vs-randomness/) | 电路下界（Shannon, Håstad, Razborov）、硬度放大、最坏情况到平均情况归约、Impagliazzo-Wigderson 定理 | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-list-decodable-codes](mini-list-decodable-codes/) | 列表可解码码、Johnson/Griesmer 界、Goldreich-Levin 定理、Reed-Solomon/Reed-Muller 码、折叠 Reed-Solomon 码 | MIT 6.841, Stanford CS254, CMU 15-855 |
| [mini-nisan-wigderson-prg](mini-nisan-wigderson-prg/) | Nisan-Wigderson PRG 构造、组合设计、混合论证、Yao XOR 引理、硬度放大、BPP 去随机化 | MIT 6.841, Stanford CS254, Berkeley CS278 |
| [mini-pseudorandom-generators](mini-pseudorandom-generators/) | PRG 构造、下一比特不可预测性、NIST 统计测试、硬核谓词、PRG 复合、Yao 定理 | MIT 6.841, Stanford CS254, CMU 15-855, Princeton COS522 |
| [mini-randomness-extractors](mini-randomness-extractors/) | 随机性提取器、最小熵、Leftover Hash Lemma、Trevisan 提取器、分散器、双源提取器、熵度量 | MIT 6.841, Stanford CS254, Berkeley CS278, CMU 15-855 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`docs/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明和知识图谱
- **去随机化管线** — 每个模块连接到更广泛的硬度与随机性框架（Arora-Barak 第 19–20 章）

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-bpp-rp-zpp-classes
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-randomized-derandomization/
├── mini-bpp-rp-zpp-classes/       # 随机化复杂度类：BPP、RP、ZPP、co-RP
├── mini-derandomizing-space/      # 空间界限去随机化：Savitch、NL=coNL、Reingold
├── mini-expander-graphs-construction/  # 谱扩展图、zigzag 积、Ramanujan 图
├── mini-hardness-vs-randomness/   # 电路下界、最坏到平均归约、IW 定理
├── mini-list-decodable-codes/     # 列表可解码码：Johnson 界、RS/RM 码、Goldreich-Levin
├── mini-nisan-wigderson-prg/      # NW PRG：组合设计、混合论证、XOR 引理
├── mini-pseudorandom-generators/  # PRG 基础：NIST 测试、硬核谓词、PRG 复合
└── mini-randomness-extractors/    # 提取器：最小熵、LHL、Trevisan、分散器
```

## 许可证

MIT
