---
title: A Matrix Expander Chernoff Bound 
date: 2023-02-09
speaker: Tushant Mittal
---

The talk is based on the paper, *A Matrix Expander Chernoff Bound* by Ankit Garg, Yin Tat Lee, Zhao Song, and, Nikhil Srivastava. [arXiv link](https://arxiv.org/abs/1704.03864).

## Setup
The most basic setup is as follows. We know using usual Chernoff that if $\\{X_i\\}$ are *independent* Radmacher $\pm 1$ random variables, then $ X = \sum_i X_i$ concentrates around zero with an exponentially decaying tail. Now, we wish to extend this to the setting when the variables are dependent. One way to do this is to reinterpret the above as a random walk. Let $G$ be a complete graph (with self loops) on $n$ vertices and $f: V\to \\{\pm 1\\}$ such that half of the vertices are labelled $1$ and the other are labelled $-1$. Now, sampling  an independent tuple $(X_1,\cdots, X_k)$ is equivalent to sampling a vertex, $v_1$, performing a walk for $k-1$ steps, and then outputting $(f(v_1),\cdots, f(v_k))$. Since the graph is complete, these samples are independent. However, as the graph becomes sparse, they are no longer independent. In the extreme regime when the graph just has self loops, it is just a tuple of $k$ copies of a single Rademacher variable. This reinterpretation, therefore, lets us smoothly interpolate between independence and complete dependence. What properties does the graph need such that such a Chernoff bound holds, *regardless of the labelling*?

One observation is that the graph needs to *mix well*. This can be observed by partitioning the vertices into two sets $V = S\cup T$ of equal size. Label $f(v) = 1$ if $v \in S$ and $-1$ if $v\in T$. This labelling has mean zero, but intuitively the walk needs to mix these two sets well for a Chernoff-like concentration to hold. Naturally then, one asks if concentration holds for an *expander graph*. And the answer, quite amazingly, is yes.  

### Formal Definitions

Let $ G = (V,E)$ be a fixed $ d$-regular graph with a normalized adjacency matrix $\rA_G$. Let the eigenvalues of $\rA_G$ be $ \lambda_1\geq \lambda_2\cdots\geq \lambda_n$. We define $\lambda(G) = \max \\{ \lambda_2, |\lambda_n| \\}$. $ G$ is an expander if $\lambda < 1$. Let $ f: V \to \mathcal{H}_d(\C)$ map vertices to $ d\times d$ Hermitian matrices such that for each $v$, $\norm{f(v)} \leq 1$ and over the uniform distribution over vertices, $\Ex{v}{f(v)} = 0$. Here, and throughout the talk, $\norm{\rA}$ denotes the spectral norm, i.e., the largest singular value. 

Let $\cW_k$ be the collection of $(k-1)$-length walks $(v_1,\cdots, v_k)$, on $G$. We will refer to it also as a distribution which is uniform over the tuples. Another way to generate this tuple is to pick $v_1$ uniformly at random and then use the transition matrix, $ A_G$, to perform a walk. 

### Statements
Let $G$ be an expander graph such that $\lambda(G)\leq \lambda$ and $\cW_k$ be the uniform distrution of tuples of $(k-1)$-length walks on $G$. Then for any $ f: V\to H_d(\C)$
 such that for each $v$, $\norm{f(v)} \leq 1$ and over the uniform distribution over vertices, $\Ex{v}{f(v)} = 0$. Then, 
\\[
\begin{align}
\small{[\text{Scalar Expander Chernoff}]}\\;\\;\\; \Pr_{(v_1,\cdots, v_k)\sim \cW_k} \\left[\abs{\sum_{i=1}^k f(v_i)} \geq \epsilon k \\right] ~\leq~ 2\mathrm{exp}\\left(\Omega(-k^2\epsilon)\\right)  \\\
\small{[\text{Matrix Expander Chernoff}]}\\;\\;\\; \Pr_{(v_1,\cdots, v_k)\sim \cW_k}  \\left[\norm{\sum_{i=1}^k f(v_i)} \geq \epsilon k \\right] ~\leq~ 2\mathrm{exp}\\left(\Omega(-k^2\epsilon)\\right)  \\\
\end{align}
\\]

## Outline
Recall that in the proof of regular Chernoff, we analyze the moment generating function, $\ex{e^{\theta \sum_i X_i}}$. There are essentially two main things, we need to handle. 
1. **Break the sum into a product**. This is an issue only in the matrix case as, multiplicativity of the exponential fails ($\small{e^{A+B} \neq e^A\cdot e^B}$). As before, this will be handled using a new trace inequality which is a multi-matrix generalization of the Golden-Thompson inequality.

2. **Dependence of the random Variables.**  This is a novel constraint that we haven't encountered so far. We handle this first via Healy's method which rewrites the expectation as a linear map.

<!---
### Template Chernoff

Let $X = \sum_i X_i$ where $X_i$ are mean zero random variables. Fix $t > 0\in \R$.
\\[
\begin{align}
\prob{\tr(X) > \epsilon k} ~&=~ \prob{e^{t X} > e^{\epsilon kt}}\\;\\;\small{\text{ (Monotonicity of exponential)} }\\\
 &\leq e^{- \epsilon kt}\\,\ex{ \tr(e^{t X})} \\;\\;\small{\text{ (Markov)} }\\\
 ~&\leq~ e^{-\epsilon kt}\\,\ex{ \tr \prod_i  g(t X_i) } \\;\\; \small{\text{(Trace Inequality)}} \\\
 ~&\leq~ e^{-\epsilon kt}\\,\norm{(\rM_g \rA_G)^k}_{\mathrm{op}} \\;\\; \small{\text{ (Healy's)} }
\end{align}
\\]

As usual, the last step is to optimize for $t$.
--->

## Healy's Method to handle dependence

 We typically use idependence to break the term $\ex{ \prod_i e^{\theta X_i}} = \prod_i \ex {e^{\theta X_i}}$. This trick we no longer have access to. So instead, we use linear algebra to bound the term $\Ex{(v_1,\cdots, v_k)\sim \cW_k}{\prod_i{e^{\theta f(v_i)}}}$. We will write it in a more general way to enable reusing it later. 
 
 Let $g: V \to \mathcal{U}(\mathcal{H})$ be a labelling of the vertices with unitary operators acting on an inner product space $\mathcal{H}$. We will work with the vector space $U =  \C[V]\otimes \mathcal{H} $ where $\C[V]$ is a vector space with the vertices forming a formal basis. The space $U$ inherits a natural inner product from each of the constituent spaces.  Define $A_G' = A_G\otimes I$ and let $\ones$ be the normalized all-ones vector in $\C[V]$. Clearly then, $\C[V] = V^\perp\oplus V^\parallel$ where $V^{\parallel} = \span (\ones)$. This induces the following decomposition, $U = U^\parallel\oplus U^\perp$ where $U^\parallel = V^\parallel \otimes \mathcal{H}$ which corresponds to the $1$-eigenspace of $A_G'$. Define the linear operator $M_g$, 

\\[\rM_g: U\to U,\\;\\; (v\otimes x) \mapsto v\otimes (g(v)\cdot x ). \\] 

This is a clearly a linear operator. We will usually take $\mathcal{H}$ to be a space of matrices and the map will be matrix multiplication.  Let $\rP$ be a projection map to the second space, i.e., $\rP(v\otimes x) \mapsto x$. Now, the claim is, 
<!--is a linear map. (To be precise, one should write $v\otimes T_{g(v)} x$ where $T_A$ is the linear map on $\mathcal{H}$ that acts by multiplication on the right by $A$.)-->
\\[ \Ex{(v_1,\cdots, v_k)\sim \cW_k}{g(v_1)\cdots g(v_k)}
  =   \rP (\rM_gA_G')^k (\ones\otimes \rI) \\]

Formally, this can be easily proved via induction. The key idea is that the random walk is performed on the first space and as the walk progresses, the second space is used to record the product of terms obtained so far. Say we started from a vertex $v_1$ initially with an identity matrix $\rI$ . Then,
\\[ v_1\otimes I \xrightarrow{\rM_G} v_1\otimes g(v_1)  \xrightarrow{\rA_G'} \left(\Ex{v_2\sim v_1}{v_2}\right) \otimes g(v_1) \xrightarrow{\rM_G} \Ex{v_2\sim v_1}{v_2 \otimes g(v_1)g(v_2)}
    \\]

Therefore, 
    \\[ \norm{\Ex{(v_1,\cdots, v_k)\sim \cW_k}{g(v_1)\cdots g(v_k)}} \leq  \norm{\rP(\rM_g\rA_G')^k(\ones\otimes \rI)}  \\]


### The Trace Inequality 
\\[
\begin{align}
\tr(e^{t X}) ~&\leq~ d^{1-\pi/4} \int_{-\pi/2}^{\pi/2}  \tr\\, \left( \prod_{j=1}^k \mathrm{exp}\left(\frac{4te^{i\phi}}{\pi}X_j\right) \prod_{j=k}^1 \mathrm{exp}\left(\frac{4te^{-i\phi}}{\pi}X_j  \right)    \right) d\mu(\phi)\\\
\ex{\tr(e^{t X})} ~&\leq~ d^{1-\pi/4} \int_{-\pi/2}^{\pi/2}  \tr\\, \Ex{(v_1,\cdots, v_k)\sim \cW_k}{\left( \prod_{j=1}^k \mathrm{exp}\left(\frac{4te^{i\phi}}{\pi}X_j\right) \prod_{j=k}^1 \mathrm{exp}\left(\frac{4te^{-i\phi}}{\pi}X_j  \right)}    \right) d\mu(\phi
\end{align}
\\] 

### Finishing the proof

Based on the trace inequality, we set $\mathcal{H} = \C^{d\times d}$ with the Hilber-Schmidt inner product $\ip{A}{B} = \tr (AB^*)$ and $g$ is defined as the following linear operator acting on left and right,
\\[
    g_{\phi,t}(v)\cdot \rM = \mathrm{exp}\left(\frac{4te^{i\phi}}{\pi} f(v) \right)\cdot  \rM \cdot \mathrm{exp}\left(\frac{4te^{-i\phi}}{\pi}f(v)  \right) 
\\]

Therefore, the term inside the integral then is $\ip{\rI_d}{\Ex{(v_1,\cdots, v_k)\sim \cW_k}{g(v_1)\cdots g(v_k)} \rI_d}$. Therefore, it suffices to bound $\opnorm{(\rM_g\rA_G')^k}\leq\opnorm{\rM_g\rA_G'}^k $. This operator norm computation is carried out in Lemma 4.4, Claim 4.5, and, Claim 4.6 in the paper. We write out the main result,

{{< css.inline >}}
<span class="theorem">(Norm bounds) <i> For any $\phi,t$, $\norm{\rP(\rM_g\rA_G')^k (\ones\otimes \rI)} \leq e^{ kt^2\left(1+\frac{8}{1-\lambda}\right)}$.   </i>
</span>
{{</ css.inline >}}



## "Deriving" the Trace Inequality — Multi-Matrix Golden Thompson

A first version of this inequality was first obtained in *Multivariate Trace Inequalities*
by David Sutter, Mario Berta, and Marco Tomamichel, [arxiv](https://arxiv.org/abs/1604.03023) . Therein, the RHS was computed as an integral over all of $\R$. This paper converts that to an integral over the half-circle.

We will outline a sketch of the key steps that go into the proof. The interested reader is encouraged to consult the paper which is very well-written.

1. **Recasting as a norm** — Rewrite $\log \tr e^{\sum_i X_i}$ as the Schatten-norm of a product. More precisely,
 \\[\log \tr e^{\sum_i X_i} = \lim_{\theta\to 0} \frac{2}{\theta}\log \lVert G(\theta)\rVert_{\frac{2}{\theta}}
 \\] where $G(\theta) = \prod_{i=1}^k e^{\theta X_i} $.
  
1. **Complex Interpolation Theory** — To estimate the above norm, they use ideas from interpolation theory. The basic idea is to something like Riesz-Thorin style of inequality that in it's simplest form, is the familiar Holder's inequality.
It says that if for a (bounded) linear operator $T$, one knows $\lVert T\rVert_{p_0}, \lVert T\rVert_{{p_1}}$, then one can bound all intermediate norms. In particular, if $p_0 = \infty, p_1 = 2$, then by Riesz-Thorin, $\log \norm{T}_{2/\theta} \leq \frac{1}{\theta}\log\norm{T}_2$. This is very similar to the expression we want, but our map $G(\theta)$ also varies with $\theta$!  Elias Stein, however, had proved in his thesis that one can infact extend the proof of Riesz-Thorin to such a setup. Using Stein's interpolation directly yields something like, 

\\[ \frac{1}{\theta} \log \lVert G(\theta)\rVert_{\frac{2}{\theta}} \leq \int_{-\infty}^\infty \log 
\norm{G(1+it)} \beta_{\theta}(t) dt
    \\]


3. **Poisson Integral Formula** — To convert this to integral over a disk, the authors open up the proof of Stein-Hirschman's interpolation and reprove it on a different domain. The key tool is the following,
   
{{< css.inline >}}
<span class="theorem">(Poisson Integral Formula) <i> Let $f$ be a subharmonic function defined on the unit disk. In particular, one can take $f = \log F$ for some analytic F. Then,  </i> 

\\\[
    \log F(z)\leq \frac{1}{2\pi}\int_{\pi}^{-\pi} \log F(e^{i\phi})\frac{1-|z|^2}{|e^{i\phi}-z|^2}  d\phi
    \\\]

</span>
{{</ css.inline >}}
The remaining task is to find an $F$ such that $F$ is analytic, $|F(it)|\leq 1$, $|F(\theta)| = \lVert G(\theta)\rVert_{2} $. They construct such an $F$ following Sutter--Berta--Tomamichel. 
 

## Resources 
- Garg, Ankit, Lee, Yin T., Song, Zhao, and Nikhil Srivastava. "A Matrix Expander Chernoff Bound." [arXiv](https://doi.org/10.48550/arXiv.1704.03864).

- Sutter,David, Berta, Mario, and Marco Tomamichel,  "Multivariate Trace Inequalities."
 [arxiv](https://doi.org/10.48550/arXiv.1604.03023)

- Talk by Ankit Garg at IAS [[Youtube Link]](https://www.youtube.com/watch?v=XZN6k_NfrrA).