    using PyCall
function evaluate_model_transmission(mdl,hyperparameters,filename)
    sys=pyimport("sys")
    os=pyimport("os")
    sys.path=cat(sys.path,ENV["LUMERICAL_PYTHON_PATH"],dims=1)
    lu=pyimport("lumapi")

    hyperparameters=merge(default_hyperparameters,hyperparameters)
    fdtd=lu.FDTD(filename=filename)
    Ts=zeros(2,hyperparameters[:fpoints])
    nmonitors=length(mdl.layers)+1
    wls=fdtd.getresult("T1","T")["lambda"]
    if hyperparameters[:get_intermediate]
        Ts=zeros(nmonitors,hyperparameters[:fpoints])
        for i=1:nmonitors
            Ts[i,:]=fdtd.getresult(string("T",i-1),"T")["T"]
        end
    else
        Ts[1,:]=fdtd.getresult(string("T",0),"T")["T"]
        Ts[2,:]=fdtd.getresult(string("T",nmonitors-1),"T")["T"]
    end
    fdtd.close()
    return Ts,wls

end





