function f = phaseRetrieval(dm, uabs, PRopt)
    switch class(PRopt)
        case 'ERoptions'
            f = useER(dm, uabs, PRopt);
        case 'HIOoptions'
            f = useHIO(dm, uabs, PRopt);
        case 'MDoptions'
            f = useMD(dm, uabs, PRopt);
        otherwise
            error('Unknown phase retrieval method');
    end
end

function f = useER(dm, uabs, opt)
    solver = ER(dm, opt);
    f = solver.traceBack(uabs);
end

function f = useHIO(dm, uabs, opt)
    solver = HIO(dm, opt);
    f = solver.traceBack(uabs);
end

function f = useMD(dm, uabs, opt)
    solver = MD(dm, opt);
    f = solver.traceBack(uabs);
end