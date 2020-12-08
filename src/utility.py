### This file contains utility functions used in quiver calculations.

# Linear Algebra utility function

def proj(u,v):
    # project u onto v
    return np.inner(v,u)/np.inner(u,u) * u

def normalize(v):
    norm=np.linalg.norm(v)
    if norm < 1.0E-5: 
       raise ValueError
    return v/norm

def test_orthogonality(vectors):
    mat = np.zeros(shape=(len(vectors),len(vectors)))
    for i,u in enumerate(vectors):
        for j,v in enumerate(vectors):
            inner = np.inner(u,v)
            if inner != 0.0 or inner != 1.0:
                mat[i][j] = inner
    print("Orthogonality:")
    for l in mat:
        print(l)

def schmidt(seed_vectors, rest_vectors, dimension):
    vectors = list(seed_vectors)
    test_vectors = list(rest_vectors)
    while len(vectors) < dimension:
        try:
            test_vector = test_vectors.pop()
            for v in vectors:
                test_vector -= proj(v, test_vector)
            try:
                vectors.append(normalize(test_vector))
            except ValueError:
                pass

        except IndexError:
            raise ValueError("Could not fill the appropriate dimensional space")            
    if len(vectors) < dimension:
        raise ValueError("Could not fill the appropriate dimensional space")
    else:
        return vectors
