import re

def parse_si_bp(option):
    SI_STEPS = {
        'K': 1000,                  # strictly speaking should be k...
        'M': 1000000,
        'G': 1000000000,
        'T': 1000000000000,
    }

    option = str(option).upper().strip().replace(' ', '').replace('BP', '')

    # I'd rather not use RE for this...
    try:
        bases = re.findall('-?\d+', option)[0]
        option = option.replace(bases, '')
    except IndexError:
        raise ValueError()

    bases = int(bases)
    for char in option:
        if char in SI_STEPS:
            bases *= SI_STEPS[char]
    return bases
