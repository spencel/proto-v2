

def export(items:list, fpath:str='list.tsv'):
  with open(fpath, 'w') as f:
    for item in items:
      f.write(f'{item}\n')