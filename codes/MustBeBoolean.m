function MustBeBoolean(item)

    isbool = (item == 1) | (item == 0) | (item == true) | (item == false) | (item == 'none');
    
    if ~isbool
        eid = 'function:inputError';
        msg = sprintf('Value must be boolean!\n');
        throwAsCaller(MException(eid,msg))
    end

end