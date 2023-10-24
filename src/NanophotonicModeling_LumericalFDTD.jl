module NanophotonicModeling_LumericalFDTD
    using NanophotonicModeling
    try
        println("You are using the Lumerical FDTD extension with python API path ",ENV["LUMERICAL_PYTHON_PATH"]," . To set it, enter ENV[\"LUMERICAL_PYTHON_PATH\"]=\"PATH/TO/LUMERICAL/api/python\"")
    catch y
        @error "For the Lumerical FDTD extension, you need to specify a path to the python API of Lumerical FDTD. If you try to build models without the correct path, this tool will crash. To set it, enter ENV[\"LUMERICAL_PYTHON_PATH\"]=\"PATH/TO/LUMERICAL/api/python\""
    end
    try
        println("You are using the Lumerical FDTD extension with runner path ",ENV["LUMERICAL_RUNNING_PATH"]," . To set it, enter ENV[\"LUMERICAL_RUNNING_PATH\"]=\"PATH/TO/LUMERICAL/bin/fdtd-run-local.bat\"")
    catch y
        @error "For the Lumerical FDTD extension, you need to specify a path to the Lumerical FDTD execution script. If you try to build models without the correct path, this tool will crash. To set it, enter ENV[\"LUMERICAL_RUNNING_PATH\"]=\"PATH/TO/LUMERICAL/bin/fdtd-run-local.bat\""
    end

    default_hyperparameters=Dict(:mesh_type=>"uniform",
                                 :simulation_time=>1000e-15,
                                 :dt_stability_factor=>.2,
                                 :auto_shutoff_min=>1e-7,
                                 :pml_type=>"uniaxial anisotropic PML (legacy)",
                                 :hide=>true,
                                 :x_min_bc=>"Periodic",
                                 :y_min_bc=>"Periodic",
                                 :pml_layers=>2048,
                                :dx=>20,
                                :z_add=>[50,50],
                                :polarization_angle=>0,
                                :incidence_angle=>0,
                                :azimuth_angle=>0,
                               :λmin=>500,
                              :λmax=>1500,
                             :fpoints=>1001,
                             :getfield=>false,
                             :get_intermediatd=>true,
                             :material_loading=>Dict{String,Array{<:Real,1}}())
    include("Running.jl")
    include("Evaluating.jl")
    include("Materials.jl")

    using PyCall,CSV,DataFrames


    function add_layer(fdtd,l::PlainLayer,zpos::Real,xy_span::Array{<:Real,1},material_defs::Dict,index)
        nm=string("L",index," - ",l.material)
        fdtd.addrect(material=material_defs[l.material],
                     x=0,
                     y=0,
                     z=(zpos-l.thickness/2)*1e-9,
                     x_span=xy_span[1]*1e-9,
                     y_span=xy_span[2]*1e-9,
                     z_span=l.thickness*1e-9,
                    name=nm)
        return zpos-l.thickness
    end
    function add_layer(fdtd,l::PatternedLayer,zpos::Real,xy_span::Array{<:Real,1},material_defs::Dict,index)
        nm=string("L",index," - ",l.material)
        fdtd.addrect(material=material_defs[l.material],
                     x=0,
                     y=0,
                     z=(zpos-l.thickness/2)*1e-9,
                     x_span=xy_span[1]*1e-9,
                     y_span=xy_span[2]*1e-9,
                     z_span=l.thickness*1e-9,
                    name=nm)
        for i in eachindex(l.inclusions)
            add_block(fdtd,l.inclusions[i].geometry,l.inclusions[i].material,zpos-l.thickness/2,l.thickness,material_defs,index,i)
        end
        return zpos-l.thickness
    end
    function add_block(fdtd,g::Ellipse,material::String,z::Real,z_span::Real,material_defs::Dict,lindex,index)
        nm=string("L",lindex,"I",index," - ",material)
        fdtd.addcircle(x=0.0,
                       y=0.0,
                       z=z*1e-9,
                       make_ellipsoid=true,
                       radius=g.radius[1]*1e-9,
                       radius_2=g.radius[2]*1e-9,
                       z_span=z_span*1e-9,
                       material=material_defs[material],
                      name=nm)
    end
    function add_block(fdtd,g::Rectangle,material::String,z::Real,z_span::Real,material_defs::Dict,lindex,index)
        nm=string("L",lindex,"I",index," - ",material)
        fdtd.addrect(x=0.0,
                       y=0.0,
                       z=z*1e-9,
                       x_span=g.size[1]*1e-9,
                       y_span=g.size[2]*1e-9,
                       z_span=z_span*1e-9,
                       material=material_defs[material],
                      name=nm)
    end
    function prepare_materials(fdtd,matloading)
        for m in keys(matloading)
            add_mat(fdtd,m,matloading[m])
        end
    end

                                   
    function build_model(mdl::PeriodicLayerModel,hyperparameters::Dict,filename="")
    sys=pyimport("sys")
    os=pyimport("os")
    sys.path=cat(sys.path,ENV["LUMERICAL_PYTHON_PATH"],dims=1)
    lu=pyimport("lumapi")
        hyperparameters=merge(default_hyperparameters,hyperparameters)
        fdtd=lu.FDTD(hide=hyperparameters[:hide])
        prepare_materials(fdtd,hyperparameters[:material_loading])
        zsize=sum([l.thickness for l in mdl.layers])
        fdtd.addfdtd(dimension=("3D"),
                     x=0,
                     y=0,
                     z=-zsize/2*1e-9,
                     x_span=mdl.xy_size[1]*1e-9,
                     y_span=mdl.xy_size[2]*1e-9,
                     z_span=(zsize+hyperparameters[:z_add][1]+hyperparameters[:z_add][2])*1e-9,
                     mesh_type=hyperparameters[:mesh_type],
                     dt_stability_factor=hyperparameters[:dt_stability_factor],
                     simulation_time=hyperparameters[:simulation_time],
                     auto_shutoff_min=hyperparameters[:auto_shutoff_min],
                     pml_type=hyperparameters[:pml_type],
                     dx=hyperparameters[:dx]*1e-9,
                     dy=hyperparameters[:dx]*1e-9,
                     dz=hyperparameters[:dx]*1e-9,
                     x_min_bc=hyperparameters[:x_min_bc],
                     y_min_bc=hyperparameters[:y_min_bc],
                     pml_layers=hyperparameters[:pml_layers])
        fdtd.addplane(x=0,
                      y=0,
                      z=.5*hyperparameters[:z_add][1]*1e-9,
                      polarization_angle=hyperparameters[:polarization_angle],
                      angle_theta=hyperparameters[:incidence_angle],
                      angle_phi=hyperparameters[:azimuth_angle],
                      x_span=1.1*mdl.xy_size[1]*1e-9,
                      y_span=1.1*mdl.xy_size[2]*1e-9,
                      direction="Backward",
                      override_global_source_settings=false)
        fdtd.setglobalsource("wavelength start",hyperparameters[:λmin]*1e-9)
        fdtd.setglobalsource("wavelength stop",hyperparameters[:λmax]*1e-9)
        fdtd.setglobalmonitor("frequency points",hyperparameters[:fpoints])
        fdtd.setglobalmonitor("use source limits",true)
        zpos=0.0
        fdtd.addpower(name="T0",
                      monitor_type="2D Z-normal",
                      x=0,
                      y=0,
                      z=0,
                      output_Ex=hyperparameters[:getfield],
                      output_Ey=hyperparameters[:getfield],
                      output_Ez=hyperparameters[:getfield],
                      output_Hx=hyperparameters[:getfield],
                      output_Hy=hyperparameters[:getfield],
                      output_Hz=hyperparameters[:getfield],
                      x_span=1.1e-9*mdl.xy_size[1],
                      y_span=1.1e-9*mdl.xy_size[2])
        for l in eachindex(mdl.layers)
            zpos=add_layer(fdtd,mdl.layers[l],zpos,1.1*mdl.xy_size,hyperparameters[:material_defs],l)
            if (hyperparameters[:get_intermediate]) || (l==length(mdl.layers))
                fdtd.addpower(name=string("T",l),
                              monitor_type="2D Z-normal",
                              x=0,
                              y=0,
                              z=zpos*1e-9,
                              output_Ex=hyperparameters[:getfield],
                              output_Ey=hyperparameters[:getfield],
                              output_Ez=hyperparameters[:getfield],
                              output_Hx=hyperparameters[:getfield],
                              output_Hy=hyperparameters[:getfield],
                              output_Hz=hyperparameters[:getfield],
                              x_span=1.1e-9*mdl.xy_size[1],
                              y_span=1.1e-9*mdl.xy_size[2])
              end
        end
        fdtd.addrect(x=0,
                     y=0,
                     z=(zpos-hyperparameters[:z_add][2])*1e-9,
                     x_span=1.1*mdl.xy_size[1]*1e-9,
                    y_span=1.1*mdl.xy_size[2]*1e-9,
                    z_span=2*hyperparameters[:z_add][2]*1e-9,
                    material=hyperparameters[:material_defs][mdl.sub],
                   name="substrate")
        fdtd.addrect(x=0,
                     y=0,
                     z=(hyperparameters[:z_add][1])*1e-9,
                     x_span=1.1*mdl.xy_size[1]*1e-9,
                    y_span=1.1*mdl.xy_size[2]*1e-9,
                    z_span=2*hyperparameters[:z_add][1]*1e-9,
                    material=hyperparameters[:material_defs][mdl.sup],
                   name="superstrate")
        fdtd.save(filename)
        fdtd.close()
    end
end


