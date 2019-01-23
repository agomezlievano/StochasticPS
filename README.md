# StochasticPS
Dimensionality reduction techniques to understand the process of economic development.

These are a collection of programs (some in R, some in Python) where I explore the idea that 
economic development is a constrained process of collective learning. This implies, on the 
one hand, that economies can be represented in low dimensions, and on the other, that the 
evolution itself of economies in time follows archetypical trajectories.

The question are, thus, which type of dimensionality reduction technique we should use? How 
many archetypical trajectories are there? What are the causal agents that nudge economies 
along the developmental trajectories?

I believe the literature shows evidence that economic development follows a very specific 
process that goes under the name of ''consensus dynamics''. This makes sense, if we 
take the anthropological perspective from the field of Cultural Evolution, in which 
humans excel at imitating the behavior of prestigious agents. By this 'imitation instinct' 
knowledge can be transmitted across generations accurately, but crucially, it can accumulate. 
Hence, imitation entails collective learning.

Mathematically, consensus dynamics is modeled by the equation
$$
x_{t+1} = x_{t}P,
$$
where $x$ is a vector of characteristics across nodes (e.g., cities, countries, firms) and 
$P$ is a left-stochastic matrix. That is, a matrix that propagates a vector of probabilities 
in a Markov Chain into the future such that $p_{t+1}=P p_{t}$. 

Consensus dynamics necessarily implies that the information encoded in $x$ will cluster, because 
similar nodes will resemble each other first (''they reach consensus''). In the long term all 
nodes will have the same average value (assuming the nodes all belong to the same connected 
component).

If consensus dynamics indeed describes the process of collective learning, and by implication, 
of econonomic development, then the low-dimensional representation of economies is through 
the use of the left-eigenvectors of the similarities between the nodes. 

Since economies (such as cities or countries) are multidimensional, we must represent the system 
not through vectors but matrices. If $X_t$ is the matrix of economies times economic activities, 
a possible equation of movement is thus 
$$
X_{t+1} = X_{t} P_{t} + C_{t} X_{t},
$$
where $P_t$ are the (stochastic matrix) similarities between economic activities and $C_t$ are 
the (stochastic matrix) similarities between economies.
