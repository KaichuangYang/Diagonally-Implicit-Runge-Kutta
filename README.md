# Diagonally-Implicit-Runge-Kutta
## 1 Stiff ODEs and transients

The ability to solve an Initial Value Problems (IVP) for y(x), of the form

$$
{d(y)\over dx }=f(x,y)\quad ,\quad  y(0)=0 \tag{1}
$$

depends on the characteristic lengthscales involved in the problem. When the length- scales differ by many orders of magnitude the system is called stiff and can be very difficult to solve numerically. In such stiff problems tiny errors in the large scale be- haviour can swamp the correct transient, small scale behaviour, or vice versa. These difficulties can usually avoided by using implicit numerical solvers.

An example of such a stiff system is given by the ODE system

$$
{d\over dx}
\begin{pmatrix} 
y_1 \\ 
y_2 
\end{pmatrix} = 
\begin{pmatrix} 
-1000y_1\\
1000y_1-y_2 
\end{pmatrix} 
$$

$$
y(0)=
\begin{pmatrix}
1\\
0
\end{pmatrix}
\tag{2}$$

for which the solution is

$$
\begin{pmatrix} 
y_1 \\
y_2 
\end{pmatrix}=
\begin{pmatrix} 
e^{-1000x} \\ 
{1000\over999}(e^{-x}-e^{-1000x})
\end{pmatrix}\tag{3}
$$

As shown by **figure 1**, initially y2 jumps very rapidly from 0 to 1 before decaying very slowly. When such a stiff ODE is solved using a standard explicit algorithm it is often found that a very small step-length must be take in order to ensure stability of the result. In other words if too large step length is taken the result of the numerical method, based on the explicit algorithm, does not even come close to approximating the solution (because it is unstable). However the very small step length required by the explicit solver is a waste of computational time and so an implicit solver, which does not have these stability issues, is more appropriate.

Here we restrict our attention to equations of the form

$$
{d(y)\over dx }=Ay+b(x)\quad ,\quad  y(0)=0 \tag{4}
$$

The example above in equation (2), is of the form

$$
A=
{\begin{pmatrix} 
-a_1&0\\ 
a_1&-a_2
\end{pmatrix}} {\qquad} b=0 \tag{5}
$$

With the initial data

$$
y_0={\begin{pmatrix}
1\\ 
0
\end{pmatrix}}\tag{6}
$$

equation(5) has the solution

$$
{\begin{pmatrix} 
y_1\\ 
y_2
\end{pmatrix}}=
{\begin{pmatrix} 
e^{-a_1x}\\ 
{a_1\over a_1-a_2}(e^{-a_2x}-e^{-a_1x})
\end{pmatrix}}\tag{7}
$$

Stiff behaviour, characterised by an initial rapid transient followed by a slow decay, occurs for parameter values satisfying $a_1>>a_2>0$.

![image](https://github.com/KaichuangYang/Diagonally-Implicit-Runge-Kutta/blob/main/Example/Plot/output.png)

## 2 Algorithms

### 2.1 Explicit algorithm
We are to implement the standard explict third order Runge-Kutta method, written RK3.  Assume an evenly spaced grid with spacing h (so that $x_{n+1} = x_{n} + h$). Then, given that numerical approximation to $y(x_n)$ is $y_n$, the RK3 algorithm used to compute $y_{n+1}$ (i.e the numerical approximation to $y_{n+1}$ ) from the ODE (4) is

$$
{{y^{(1)} =y_n+h[Ay_n+b(x_n)] }, \tag{8a}}
$$

$$
{y^{(2)} ={3\over 4}y_n+{1\over 4}y^{(1)}+{1\over4}h[Ay^{(1)}+b(x_n+h)]\tag{8b}} 
$$

$$
{y_{n+1} ={1\over 3}y_n+{2\over 3}y^{(2)}+{2\over3}h[Ay^{(2)}+b(x_n+h)]}\tag{8c}
$$

where A is the matrix in (4).

### 2.2 Implicit algoritm
We are also to implement the optimal two stage third order accurate Diagonally Im- plicit Runge-Kutta method, written DIRK3. Once again assume an evenly spaced grid with spacing h  h (so that $x_{n+1} = x_{n} + h$). Then, given that numerical approximation to $y(x_n)$ is $y_n$, the DIRK3 algorithm used to compute $y_{n+1}$ (i.e the numerical approximation to $y_{n+1}$ ) from the ODE (4) is

$$
{[I-h\mu A]{y^{(1)} =y_n+h\mu b(x_n+h\mu) }, \tag{9a}}
$$

$$
{[I-h\mu A]y^{(2)} =y^{(1)}+h\nu [Ay^{(1)}+b(x_n+h\mu)]+ h\mu b(x_n+h\nu +2h\mu),\tag{9b}}
$$

$$
{y_{n+1} =(1-\lambda)y_n+\lambda y^{(2)}+h\gamma[Ay^{(2)}+b(x_n+h\nu +2h\mu)],\tag{9c}}
$$

where $I$ is the identity matrix, A is the matrix in (4) and the coefficients µ, ν , γ and λ are defined as

$$
\mu={1\over2}(1-{1\over 3})\tag{10a}
$$

$$
\nu ={1\over2}(\sqrt 3 -1)\tag{10b}
$$

$$
\gamma ={3\over2(3+\sqrt 3)}\tag{10c}
$$

$$
\lambda ={3(1+\sqrt 3)\over2(3+\sqrt 3)}\tag{10d}
$$

## 3 Example
### 3.1 Example 1 Moderately Stiff
The same example as Figure 1.
### 3.2 Example 2 Stiff

$$
A={
\begin{pmatrix}
-1&0&0\\
-99&-100&0\\
-10098&9900&-10000
\end{pmatrix}}  \tag{11}
$$

$$
b={\begin{pmatrix}
cos(10x)-10sin(10x)\\
199cos(10x)-10sin(10x)\\
208cos(10x)+10000sin(10x)
\end{pmatrix}}\tag{12}
$$

$$
y_0={
\begin{pmatrix} 
0\\
1\\
0
\end{pmatrix}}\tag{13}
$$

solution

$$
y={
\begin{pmatrix} 
cos(10x)-e^{-x}\\
cos(10x)+e^{-x}-e^{-100x}\\
sin(10x)+2e^{-x}-e^{-100x}e^{-10000x}
\end{pmatrix}}\tag{12}
$$
