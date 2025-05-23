# Marker-free Alignment Toolkit


## 简介

是一个使用MPI和CUDA进行 cryo-electron tomography 对齐的高效软件工具。

## 目录

- [Cryo-ET Tomoalign Toolkit](#cryo-et-tomoalign-toolkit)
  - [简介](#简介)
  - [目录](#目录)
  - [编译](#编译)
  - [依赖](#依赖)
  - [使用方法](#使用方法)
  - [参数和配置](#参数和配置)
  - [示例](#示例)



## 编译
```bash
cd build
cmake ..
make -j 8
```

### 依赖

- MPI (在aretomo/src/mrc/mrcstack.h中的mpi库引用需要改一下)
- CUDA


## 使用方法

该软件提供了多种命令行选项以便于用户能更灵活地运行 cryo-ET 对齐。基本的命令行结构如下：

```bash
./cuda [options]
```

## 参数和配置
下面是可用选项（options）的列表:

`-INPUT(-i) 输入文件名`

用于重建的 MRC 文件。


`-OUTPUT(-o) 输出文件名`

结果的 MRC 文件名。


`-TILEFILE(-t) 倾斜角度文件名`

倾斜角度。


`-GEOMETRY(-g) 六个整数`

几何信息：偏移量、倾斜轴角度、z轴偏移、厚度、投影匹配重建厚度、输出图像的缩小比例、使用的GPU（默认为0，如果只有一个GPU不用输入）。


`-AXIS 重建轴`

重建轴：y（默认）或 z。


`-METHOD(-m) 方法名称(I,R)`（还没写）

- 反投影：BPT
- 滤波反投影：FBP
- 加权反投影：WBP
- SART：SART,迭代次数,松弛参数
- SIRT：SIRT,迭代次数,松弛参数


`-help(-h)`

帮助信息。


## 示例
下面的例子展示了如何使用投影匹配厚度300来对齐：

```bash
./cuda -i /home/liuzh/tomoalign/data/BBb.st -o /home/liuzh/tomoalign/data/alignresult.mrc -a /home/liuzh/tomoalign/data/BBb.rawtlt -g 0,0,0,0,300,2,0 -m SIRT,10,0.2
```
300表示投影匹配时反投影使用的厚度，2表示输出的结果长宽各缩小两倍，最后一个0表示使用0号GPU

同时最后会输出得到的对齐参数和使用时间，可以拿来参考


