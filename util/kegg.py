def KeggIdFromInt(id):
    """Makes a string KEGG id from an integral one.
    
    Args:
        id: the integer KEGG ID.
        
    Returns:
        The string format of the KEGG ID (properly 0-padded).
    """
    return 'C%05d' % id 
