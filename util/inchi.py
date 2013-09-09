def AchiralInchi(inchi):
    """Removes chirality information from an InChI identifier.
    
    Args:
        inchi: the input InChi potentially containing chirality info.
    
    Returns:
        The identifier with chirality stripped.
    """
    fields = inchi.split('/')
    keeper = lambda f: f[0] not in ('t', 'm', 's')
    return '/'.join(filter(keeper, fields))