# 气象中的EOF和SVD详解
自己的理解，作参考之用。

### EOF的计算过程
对于一个时空数据集，观测点位于空间中m个地点——![](https://latex.codecogs.com/svg.latex?x_1,x_2,...,x_m)，观测时间是长为𝑛的时间序列——![](https://latex.codecogs.com/svg.latex?t_1,t_2,...,t_n)。这些观测值可以用一个![](https://latex.codecogs.com/svg.latex?M%5Ctimes%20N)的矩阵![](https://latex.codecogs.com/svg.latex?F)来表示，![](https://latex.codecogs.com/svg.latex?F)的行是某个地点在观测期内所有时间点的观测值，而列是某个时间点上地图上所有观测点的观测值。下面的异常矩阵![](https://latex.codecogs.com/svg.latex?A)是![](https://latex.codecogs.com/svg.latex?F)矩阵中的每个元素减去时间均值（即各行均值）而得到的。

<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?A%3D%5Cbegin%7Bbmatrix%7D%20a_%7B1%2C1%7D%20%26%20a_%7B1%2C2%7D%20%26%20%5Chdots%20%26%20a_%7B1%2Cn%7D%5C%5C%20a_%7B2%2C1%7D%20%26%20a_%7B2%2C2%7D%20%26%20%5Chdots%20%26%20a_%7B2%2Cn%7D%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cvdots%20%5C%5C%20a_%7Bn%2C1%7D%20%26%20a_%7Bn%2C2%7D%20%26%20%5Chdots%20%26%20a_%7Bn%2Cn%7D%20%5Cend%7Bbmatrix%7D"/></div>
<br/>

一个常用的计算EOF的方法是首先计算协方差矩阵。设样本对应的多为随机变量为![image](https://latex.codecogs.com/svg.latex?X=[X_1,X_2,X_3,...,X_n]^T)，样本集合为![](https://latex.codecogs.com/svg.latex?%5C%7Bx_%5Ccdot%20j%3D%5Bx_1j%2Cx_2j%2C...%2Cx_nj%5D%5ET%7C1%5Cleq%20j%5Cleq%20m%5C%7D),![](https://latex.codecogs.com/svg.latex?m)为样本数量。与样本方差的计算相似，𝑎和𝑏两个维度样本的协方差公式为![](https://latex.codecogs.com/svg.latex?1%5Cleq%20a%5Cleq%20n%2C%201%5Cleq%20b%5Cleq%20n)，![](https://latex.codecogs.com/svg.latex?n)为样本维度。

<br/>
<div align=center><img src="https://latex.codecogs.com/svg.latex?q_%7Ba%2Cb%7D%3D%5Cfrac%7B%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7Baj%7D-%5Cbar%7Bx%7D_a%29%28x_%7Bbj%7D-%5Cbar%7Bx%7D_b%29%7D%7Bm-1%7D"/></div>
<br/>

这里分母![](https://latex.codecogs.com/svg.latex?m-1)是因为随机变量的数学期望未知，以样本均值代替，自由度减一。这样，所有样本可以表示成一个![](https://latex.codecogs.com/svg.latex?n%5Ctimes%20m)的矩阵。我们以![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D)表示样本的协方差矩阵。

<div align="center"><img src="https://latex.codecogs.com/svg.latex?\hat{\Sigma&space;}=\begin{bmatrix}&space;q_{1,1}&space;&&space;q_{1,2}&space;&&space;\hdots&space;&&space;q_{1,n}\\&space;q_{2,1}&space;&&space;q_{2,2}&space;&&space;\hdots&space;&&space;q_{2,n}\\&space;\vdots&space;&&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;q_{n,1}&space;&&space;q_{n,2}&space;&&space;\hdots&space;&&space;q_{n,n}&space;\end{bmatrix}"></div>
<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?%3D%5Cfrac%7B1%7D%7Bm-1%7D%5Cbegin%7Bbmatrix%7D%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%20%26%20%5Chdots%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%5C%5C%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%20%26%20%5Chdots%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%5C%5C%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cddots%20%26%20%5Cvdots%5C%5C%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%28x_%7B1j%7D-%5Cbar%7Bx%7D_1%29%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%28x_%7B2j%7D-%5Cbar%7Bx%7D_2%29%20%26%20%5Chdots%20%26%20%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%28x_%7Bnj%7D-%5Cbar%7Bx%7D_n%29%20%5Cend%7Bbmatrix%7D"/></div>
<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?%3D%5Cfrac%7B1%7D%7Bm-1%7D%5Csum_%7Bj%3D1%7D%5E%7Bm%7D%28%5Cmathbf%7Bx%7D_%7B%5Ccdot%20j%7D-%5Cbar%7B%5Cmathbf%7Bx%7D%7D%29%28%5Cmathbf%7Bx%7D_%7B%5Ccdot%20j%7D-%5Cbar%7B%5Cmathbf%7Bx%7D%7D%29%5ET"/></div>
<br/>

接着算得特征值:
<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?R%3DP%5CLambda%20P%5E%7B-1%7D"/></div>
<br/>

### 奇异值分解（SVD）

SVD也是对矩阵进行分解，但是和特征分解不同，SVD并不要求要分解的矩阵为方阵。假设我们的矩阵![](https://latex.codecogs.com/svg.latex?A)是一个![](https://latex.codecogs.com/svg.latex?M%5Ctimes%20N)的矩阵，那么我们定义矩阵![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D)的SVD为：

<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D%3DU%5CSigma%20V%5ET"/></div>
<br/>

其中![](https://latex.codecogs.com/svg.latex?U)是一个![](https://latex.codecogs.com/svg.latex?m%5Ctimes%20m)的矩阵，![](https://latex.codecogs.com/svg.latex?%5CSigma)是一个![](https://latex.codecogs.com/svg.latex?m%5Ctimes%20n)的矩阵，除了主对角线上的元素以外全为0，主对角线上的每个元素都称为奇异值，![](https://latex.codecogs.com/svg.latex?V)是一个![](https://latex.codecogs.com/svg.latex?n%5Ctimes%20n)的矩阵。![](https://latex.codecogs.com/svg.latex?m%5Ctimes%20U)和![](https://latex.codecogs.com/svg.latex?m%5Ctimes%20V)都是酉矩阵，即满足![](https://latex.codecogs.com/svg.latex?U%5ETU%3DI%2CV%5ETV%3DI)。下图可以很形象的看出上面SVD的定义：

<div align="center"><img src="https://user-images.githubusercontent.com/76199161/138597327-5021e994-b3fb-4f60-bb50-f2dbbc0f3e3c.png"/></div>
<div align="center">图1</div>

接着我们计算![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D%5ET%5Chat%7B%5CSigma%7D)的n个特征值对应的n个特征向量v和![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D%5Chat%7B%5CSigma%7D%5ET)的m个特征值对应的m个特征向量u。由于![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%20%7D)除了对角线上是奇异值其他位置都是0，那我们只需要求出每个奇异值σ就可以了。

我们注意到:
<br/>
<div align="center"><img src="https://latex.codecogs.com/svg.latex?%5Chat%7B%5CSigma%7D%3DU%5CSigma%20V%5ET%5Crightarrow%20%5Chat%7B%5CSigma%7D%3DU%5CSigma%20V%5ETV%5Crightarrow%20%5Chat%7B%5CSigma%7DV%3DU%5CSigma%20%5Crightarrow%20%5Chat%7B%5CSigma%7Dv_i%3D%5Csigma%20_iu_i%5Crightarrow%20%5Csigma%20_i%3D%5Cfrac%7B%5Chat%7B%5CSigma%7Dv_i%7D%7Bu_i%7D"/></div>
<br/>

这样我们可以求出我们的每个奇异值，进而求出奇异值矩阵Σ。
