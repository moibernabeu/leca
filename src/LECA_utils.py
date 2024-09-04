from utils import count_list_items
from tree_utils import read_lineage

num2words = {1: 'One', 2: 'Two', 3: 'Three', 4: 'Four', 5: 'Five',
             6: 'Six', 7: 'Seven', 8: 'Eight', 9: 'Nine', 10: 'Ten'}

LECA_criteria = {'stems': 2,
                 'supergroups': [3, 5, 7],
                 'sp': 5}

lng = read_lineage('/gpfs/projects/bsc40/current/mgil/documents/leca_V3/data/toldb.lng')


def get_sp(header):
    if 'Prok|' in header:
        if 'rvdb' in header:
            return header.replace('Prok|', '').split('_', 1)[0]
        else:
            return header.split('_', 2)[1]
    else:
        return header.split('_', 1)[0]


def is_LECA(ids: list, lng: dict=lng, criteria: dict=LECA_criteria):
    sp = [get_sp(x) for x in ids]

    group_data = {}
    group_data['sp'] = len(set(sp))
    group_data['stems'] = len(set([lng[x]['stem'] for x in sp]))
    group_data['supergroups'] = len(set([lng[x]['sg'] for x in sp]))

    group_criteria = {}
    for criterion, v in criteria.items():
        if isinstance(v, list):
            for subcrit in v:
                group_criteria[f'{num2words[subcrit]}_{criterion}'] = group_data[criterion] >= subcrit
        else:
            group_criteria[f'{num2words[v]}_{criterion}'] = group_data[criterion] >= v
    
    min_sg = group_criteria[f'{num2words[min(LECA_criteria["supergroups"])]}_supergroups']
    min_sp = group_criteria[f'{num2words[LECA_criteria["sp"]]}_sp']
    min_st = group_criteria[f'{num2words[LECA_criteria["stems"]]}_stems']

    is_leca_set = min_sg and min_sp and min_st

    group_criteria = {**group_criteria, 'is_LECA': is_leca_set, **group_data}

    return group_criteria
