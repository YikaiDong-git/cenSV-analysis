"""Enumerate all call subsets with a required net edit length, within a state cap."""


class SubsetSearchLimit(RuntimeError):
    pass


def net_bounds(calls, indices):
    nets=[len(calls[i]['edit'][2])-calls[i]['edit'][1]+calls[i]['edit'][0] for i in indices]
    return sum(x for x in nets if x<0),sum(x for x in nets if x>0)


def compatible_subsets(calls, indices, target, states, max_states):
    indices=list(indices); n=len(indices)
    nets=[len(calls[i]['edit'][2])-calls[i]['edit'][1]+calls[i]['edit'][0] for i in indices]
    lower=[[float('inf')]*(n+1) for _ in range(n+1)]
    upper=[[float('-inf')]*(n+1) for _ in range(n+1)]
    lower[n][0]=upper[n][0]=0
    for pos in range(n-1,-1,-1):
        lower[pos][0]=upper[pos][0]=0
        for count in range(1,n-pos+1):
            lower[pos][count]=min(lower[pos+1][count],nets[pos]+lower[pos+1][count-1])
            upper[pos][count]=max(upper[pos+1][count],nets[pos]+upper[pos+1][count-1])

    def visit(pos,count,needed,chosen):
        states[0]+=1
        if states[0]>max_states:
            raise SubsetSearchLimit('Net-compatible subset enumeration reached its state cap')
        if count>n-pos or needed<lower[pos][count] or needed>upper[pos][count]:
            return
        if count==0:
            if needed==0:yield tuple(chosen)
            return
        yield from visit(pos+1,count-1,needed-nets[pos],chosen+[indices[pos]])
        yield from visit(pos+1,count,needed,chosen)

    for count in range(1,n+1):
        yield from visit(0,count,target,[])
