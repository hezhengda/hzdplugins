# This document contains many small functions that can help us find relevant information about each node

import pandas as pd

def uuid(pk):

    """

    Give the uuid of certain node

    Parameters:

    pk:
        The pk of the node

    Return: uuid of that node

    """

    node = load_node(pk)
    print(node.uuid)

def showresults(results):

    """

    show results in pandas form.

    Parameters:

    results:
        The results dictionary that contains all relevant computational information.

    Return: A panda object that we can use to display the dataframe 

    """

    df = pd.DataFrame.from_dict(results, orient='index')
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    return df
