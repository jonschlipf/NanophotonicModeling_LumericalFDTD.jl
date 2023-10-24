using CSV,DataFrames
function add_mat(fdtd,name,data)
    fdtd.setmaterial(fdtd.addmaterial("Sampled 3D data"),"name",name)
	fdtd.setmaterial(name,"sampled data",data)
end
