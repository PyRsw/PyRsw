# Helper for formatting time strings
def smart_time(t):
    tstr = 't = '

    if t < 2*60.:
        tstr += '{0:.4f} sec'.format(t)
    elif t < 90.*60:
        tstr += '{0:.4f} min'.format(t/60)
    elif t < 48.*60*60:
        tstr += '{0:.4f} hrs'.format(t/(60*60))
    else:
        tstr += '{0:.4f} day'.format(t/(24*60*60))

    return tstr
