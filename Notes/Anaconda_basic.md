# Anaconda使用入门

## 简介

 - `Conda`是一个开源的包、环境管理器，可以用于在同一个机器上安装不同版本的软件包及其依赖，并能够在不同的环境之间切换
 - `Anaconda`包括`Conda`、`Python`以及一大堆安装好的工具包，比如：`numpy`、`pandas`等
 - `Miniconda`包括`Conda`、`Python`  

一般我们下载使用的就是`Anaconda`，包括了基本的一些工具包，
`conda`就是用于管理包和环境的命令行工具
`pip`也是用于管理python包的工具，仅可以安装python packages

二者区别：conda能别pip做更多工作，它可以安装package依赖的其他包，而pip则不可以。conda安装的packages是a new packaging format，不能使用pip安装。conda可以看出是pip的升级。

## `conda`安装

```
# 下载 anaconda3,右键获取链接 wget + https

# sh anaconda.sh 

# 设置安装目录
```

## `conda`管理

```
# 确认conda已安装
conda --version

# 更新conda版本
conda update conda
```
其中更新命令不仅仅会更新conda的版本，同时会自动更新相关的包，
其实，我们也可以使用这个命令来更新Anaconda版本
```
conda update anaconda
```



### 修改anaconda镜像

```
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
```



## 环境管理

这里的环境指的是不同的软件版本及其依赖所构成的环境，
环境之间“绝缘”，相同软件包的不同版本可以存在于同一机器下  

```
# 创建新环境
conda create --name snowflakes biopython
```
其中snowflakes代指环境的名称，biopython指要在新环境中添加的软件包，
这里并没有指定新的环境所要使用的Python版本，所以会使用当前环境使用的Python版本

```
# 查看当前环境
conda info --envs
# conda environments:
#
# root                  *  C:\Program Files\Anaconda3
# snowflakes               C:\Program Files\Anaconda3\envs\snowflakes
```
上述命令会列出当前所有可用的环境及其路径，并在当前使用的环境前添加`*`  

root是在安装Anaconda时自动创建的环境名称，
其Python版本根据选择的Anaconda版本而定

```
# 创建环境时指定Python版本
conda create --name bunnies python=3 astroid babel
```
在创建环境指定软件包时，可以使用`package_name=version_number`
的方式来指定要使用的软件版本

```
# 切换环境
# Linux, OSX: 
# source activate snowflakes
#
# Windows:
activate snowflakes

# 切换回默认环境(root)
# Linux, OSX: 
# source deactivate
#
# Windows:
deactivate
```

其实，还可以复制一个和指定环境完全相同的环境，
只要在创建时添加`--clone`参数指定相应的环境名称即可
```
# 复制环境
conda create --name flowers --clone snowflakes
```
另外，环境也可以在不同机器之间进行复制，
只要将要复制的环境导出为`*.yml`配置文件，
再到指定机器上创建时指定配置文件即可
```
# 导出配置文件
conda env export --name snowflakes > snowflakes.yml

# 根据配置文件导入环境
conda env create -f snowflakes.yml
```

## 软件包管理

```
# 查看所有已安装的软件包
conda list
```
可用的完整软件包列表可以在<http://docs.continuum.io/anaconda/pkg-docs.html>中查找，
所有的软件包都按照Python的版本进行了分类  

当我们想要安装某个软件包时，可以直接在命令行中进行查找并安装

```
# 查找软件包
# 罗列出所有可用的版本并在已经安装的版本前加*
conda search beautifulsoup4

# 安装软件包
conda install --name beautifulsoup4=4.4.1
```

另外，也可以<http://anaconda.org>网站上搜索想要的软件包，
根据页面上的提示执行相应的命令即可安装  

最后，同样的可以使用`pip`命令来安装软件包

```
pip install XXX
```

而更新软件包可以使用`update`命令
```
conda update --name snowflakes beautifulsoup4=4.5.1
```

## `python`管理

对于`conda`来说，其实python也是一个软件包，
所以，`python`的管理基本和软件包管理相同

```
# 查找可用python版本
conda search --full-name python
```
查找名称完全匹配python的软件包，而不是名称还有python的软件包，
可以在创建环境时指定python版本
```
conda create -n snakes python=3.4
```

## 卸载包、环境

```
# 卸载包
# 删除指定环境中的指定包
conda remove --name snowflakes biopython

# 卸载环境
# --all参数表示移除环境中的所有软件包，即删除整个环境
conda remove --name snakes --all
```

> TIPS:  
> 所有命令都可以使用`--help`参数来查找详细的参数说明及用法

> 参考链接：  
> <http://conda.pydata.org/docs/test-drive.html>  
> <https://docs.continuum.io/_downloads/Anaconda_CheatSheet.pdf>
> https://github.com/search?q=baiyangcao



## pip基础用法

#### 1 pip 安装包

* pip install packagenames

#### 2 pip查看已安装的包

* pip show --files package

#### 3 pip检查哪些包需要更新

* pip list --outdated 

#### 4 pip升级包

* pip install --upgrade packages

#### 5 pip卸载包

* pip uninstall packagenames

#### 6. pip参数解释

> pip --help
> Usage:   
>   pip <command> [options]
> Commands:
>   install                     安装包.
>   uninstall                   卸载包.
>   freeze                      按着一定格式输出已安装包列表
>   list                        列出已安装包.
>   show                        显示包详细信息.
>   search                      搜索包，类似yum里的search.
>   wheel                       Build wheels from your requirements.
>   zip                         不推荐. Zip individual packages.
>   unzip                       不推荐. Unzip individual packages.
>   bundle                      不推荐. Create pybundles.
>   help                        当前帮助.
> General Options:
>   -h, --help                  显示帮助.
>   -v, --verbose               更多的输出，最多可以使用3次
>   -V, --version               现实版本信息然后退出.
>   -q, --quiet                 最少的输出.
>   --log-file <path>           覆盖的方式记录verbose错误日志，默认文件：/root/.pip/pip.log
>   --log <path>                不覆盖记录verbose输出的日志.
>   --proxy <proxy>             Specify a proxy in the form [user:passwd@]proxy.server:port.
>   --timeout <sec>             连接超时时间 (默认15秒).
>   --exists-action <action>    Default action when a path already exists: (s)witch, (i)gnore, (w)ipe, (b)ackup.
>   --cert <path>               证书