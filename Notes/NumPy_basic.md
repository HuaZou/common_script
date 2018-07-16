
## 关于NumPy

* `Python 3.5`

* Jupyter Notebook 

* NumPy 是python的矩阵类型，提供大量矩阵处理函数，内部通过C实现

* 包含两种数据类型，数组array和矩阵matrix

### 构建数组array

* 通过tuple构建array


```python
from numpy import *
tuple_elements = (4, 5, 6)
array_tuple = array(elements)
array_tuple
```


    array([4, 5, 6])



* 通过list构建array


```python
list_elements = [1, 2, 3]
array_list = array(list_elements)
array_list
```


    array([1, 2, 3])



* 构建多维array


```python
list1 = list(range(1, 4))
list2 = list(range(4, 7))
multip_array = array([list1, list2])
multip_array
```


    array([[1, 2, 3],
           [4, 5, 6]])



### array基本操作

* index索引


```python
array = multip_array
array[0][1]
```


    2



* arry的对应相乘


```python
arr2 = array*2
arr2

arr3 = array*array
arr3
```


    array([[ 1,  4,  9],
           [16, 25, 36]])



### 构建矩阵matrix


```python
# tuple构建
matrix(tuple_elements)
```


    matrix([[4, 5, 6]])




```python
# list构建
matrix(array_list)
```


    matrix([[1, 2, 3]])




```python
# array构建matrix
matrix(array)
```


    matrix([[1, 2, 3],
            [4, 5, 6]])



### matrix基本操作


```python
matrix_test = matrix(array)

# 查看维数
matrix_test.shape
```


    (2, 3)




```python
# 取值，分片方法
matrix_test[1, 2] # 第二行第三列元素
```


    6




```python
matrix_test[1, :] # 取第二行所有列
```


    matrix([[4, 5, 6]])




```python
# 转置
matrix_test.T
```


    matrix([[1, 4],
            [2, 5],
            [3, 6]])




```python
# 乘法，规则: (m, n)*(n,p)
matrix_test*(matrix_test.T)
```


    matrix([[14, 32],
            [32, 77]])




```python
# 对应元素相乘  multiply(a,b)
multiply(matrix_test, matrix_test)
```


    matrix([[ 1,  4,  9],
            [16, 25, 36]])




```python
# 排序sort，注意是原地排序，会改变原始数据，这个和pandas中的index操作不一样
mat1 = matrix([[1, 2, 4, 3]])
mat1.sort()
mat1
```


    matrix([[1, 2, 3, 4]])




```python
# 获得矩阵中每个元素的排序序号
mat2 = mat([[3,1,4],[2,3,4]])
mat2.argsort()

```


    matrix([[1, 0, 2],
            [0, 1, 2]], dtype=int64)



***The main advantage of numpy arrays is that they are more general than 2-dimensional matrices. What happens when you want a 3-dimensional array? Then you have to use an ndarray, not a matrix object. Thus, learning to use matrix objects is more work -- you have to learn matrix object operations, and ndarray operations.***
