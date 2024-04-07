
import math

import numpy


class Random():
  
  @staticmethod
  def get_balanced_randomization(
    items: list|set,
    qty: int = 10,
  ) -> list:
    """Uses: https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.choice.html

    Args:
        items (list | set): _description_
        qty (int): The size of the randomized data set.


    Returns:
        list: _description_
    """
    item_qty = len(items)

    rng = numpy.random.default_rng()
    iterations = math.ceil(qty / item_qty)

    res = list()

    # Generate random list
    for i in range(iterations):
      rng.shuffle(items)
      res.extend(items)
    
    # Get qty of items to prune
    remove_qty =  qty % item_qty
    if remove_qty > 0:
      remove_qty = item_qty - remove_qty 

    # Shuffle entire list of items for final randomization of entire sequence
    # of events
    rng.shuffle(res)

    # Remove trailing items if items/qty is not a whole number
    for i in range(remove_qty):
      res.pop()
    
    return res
  

  @staticmethod
  def _test_get_balanced_randomization(
    items: list|set,
    qty: int = 10,
  ) -> list:
    """This methoud should be deleted and put in the tests instead.

    Args:
        items (list | set): _description_
        qty (int, optional): _description_. Defaults to 10.

    Returns:
        list: _description_
    """
    
    res = Random.get_balanced_randomization(items, qty)
    print(res)

    meta = dict()
    for item in items:
      meta[item] = 0
    
    for item in res:
      meta[item] += 1
    
    return meta