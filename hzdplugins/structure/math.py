import numpy as np

def rotation_matrix_from_vectors(vec1, vec2):

    """
    Find the rotation matrix that aligns vec1 to vec2
    This function is taken from Peter's answer on stackoverflow:
    `Question Link <https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d
    -space>`_
    Parameters:
    vec1:
        Numpy array. A 3d "source" vector
    vec2:
        Numpy array. A 3d "destination" vector
    Return:
        A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """

    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)

    tmp = 0

    for item, num_a in enumerate(a):
        num_b = b[item]
        if (num_a != 0) and (num_b != 0):
            tmp += num_a / num_b

    tmp = tmp / 3 # calculate the average

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    if abs(s) < 0.0000000000001: # which means vec1 and vec2 are on the same line, so their cross is 0
        if tmp > 0:
            return [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]]
        else:
            rotation_matrix = [[-1, 0, 0],
                               [0, -1, 0],
                               [0, 0, -1]]
            return rotation_matrix

    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix

def rotation_matrix_euler(alpha, beta, gamma):

    """

    :code:`rotation_matrix_euler` can return a matrix which can rotate a vector in three directions.

    :param alpha: rotate along axis 1 (arbitrary, you can image it)
    :type alpha: radians

    :param beta: rotate along axis 2
    :type beta: radians

    :param gamma: rotate along axis 3
    :type gamma: radians

    :returns: The rotational matrix that can help to rotate a certain vector in 3D space.
    :rtype: np.array matrix

    """

    rot_mat_a = np.array([[1, 0, 0],
                          [0, np.cos(alpha), np.sin(alpha)],
                          [0, -np.sin(alpha), np.cos(alpha)]])

    rot_mat_b = np.array([[np.cos(beta), 0, -np.sin(beta)],
                          [0, 1, 0],
                          [np.sin(beta), 0, np.cos(beta)]])

    rot_mat_c = np.array([[np.cos(gamma), np.sin(gamma), 0],
                          [-np.sin(gamma), np.cos(gamma), 0],
                          [0, 0, 1]])

    tmp = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])

    tmp = np.matmul(tmp, rot_mat_a)
    tmp = np.matmul(tmp, rot_mat_b)
    tmp = np.matmul(tmp, rot_mat_c)

    return tmp