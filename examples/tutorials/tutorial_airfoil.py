"""NACA 0012 CST fit shared by the Philote tutorials.

This is the multipoint example's initial degree-4 thickness polynomial,
elevated to degree 5 without changing the shape. Coordinates use unit chord.
"""

import numpy as np

CST_0012 = np.array(
    [
        0.1706635358053721,
        0.15626369325824663,
        0.1546327069833401,
        0.14462202838023833,
        0.13318428664336585,
        0.15537328876198608,
    ]
)
TE_UPPER = 0.00063
TE_LOWER = -0.00063
