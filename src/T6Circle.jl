module T6Circle
    # Text added
    using RandomWalk
    using ContinuousTimeMarkovChain
    using ProbabilityKillZone
    using ContinuousTimeMarkovChain: State_1, State_2, State_3
    import ContinuousTimeMarkovChain as ctm
    using ProbabilityKillZone
    using RandomWalk
    import RandomWalk as rw
    import RandomWalk: time_series
    using JLD2  
    using Distributions
    using Random
    using LinearAlgebra
    using Dates
    using Colors
    using StatsBase

 
    export 
        AbstractT6SSMachineModel,
        StateDependentRandomWalkParameters,
        MultipleMachinesStateDependentRandomWalk,
        KillZone,
        CapsuleRadius,
        CapsuleOverlap,
        CapsuleKillThreshold,
        set_λ₀!,
        get_λ₀,
        get_num_machines,
        set_noise!,
        get_noise,
        set_λ₁!,
        get_λ₁,
        get_angle,
        file_name_saver,
        set_file_name_path!,
        rta,
        get_estimated_waiting_time,
        scale_with_state
        





    ######################################################
    # Kill angle struct and functions 
    ######################################################
    // # functions 
        abstract type AbstractKillZone end
        struct KillAngle <: AbstractKillZone
            angle :: Float64
        end

        mutable struct CapsuleRadius <: AbstractKillZone
            radius :: Union{Int,Float64}
        end

        mutable struct CapsuleOverlap <: AbstractKillZone
            overlap :: Float64
        end

        mutable struct CapsuleKillThreshold <: AbstractKillZone
            threshold :: Float64
        end


        mutable struct KillZone <: AbstractKillZone
            radius :: CapsuleRadius
            overlap :: CapsuleOverlap
            kill_threshold :: CapsuleKillThreshold
        end

        function get_capsule_radius(kz::KillZone) 
            kz.radius.radius
        end

        function set_capsule_radius!(kz::KillZone,radius)
            kz.radius.radius = radius
        end

        function get_capsule_overlap(kz::KillZone)
            kz.overlap.overlap
        end

        function set_capsule_overlap!(kz::KillZone,overlap)
            kz.overlap.overlap = overlap
        end

        function get_capsule_kill_threshold(kz::KillZone)
            kz.kill_threshold.threshold
        end

        function set_capsule_kill_threshold!(kz::KillZone,threshold)
            kz.kill_threshold.threshold = threshold
        end

        function get_angle(ka::KillAngle)
            ka.angle
        end


        function get_kill_angle(kz::KillZone)
            r = get_capsule_radius(kz)
            o = get_capsule_overlap(kz)
            γ = get_capsule_kill_threshold(kz)
            a = r
            b = 2*r - o
            c = γ*r
            KillAngle(acos(((a^2 + b^2) - c^2)/(2*a*b)))
        end

        function get_angle(kz::KillZone)
            get_kill_angle(kz) |> get_angle
        end
    // # end

    ######################################################
    # Struct and types
    ######################################################
    // # functions
        abstract type AbstractT6SSMachineModel end
        struct RandomWalkOnly <: AbstractT6SSMachineModel end
        struct StateTransitionProcessOnly <: AbstractT6SSMachineModel end
        struct RandomWalkOnlyPerRandomDrawλ₁ <: AbstractT6SSMachineModel end
        struct StateIndependentRandomWalk <: AbstractT6SSMachineModel end
        struct StateDependentRandomWalk <: AbstractT6SSMachineModel end
        struct MultipleMachinesRandomWalkOnly <: AbstractT6SSMachineModel end
        struct MultipleMachinesStateDependentRandomWalk <: AbstractT6SSMachineModel end



        mutable struct StateDependentRandomWalkParameters
            end_time :: Union{Int,Float64}
            time_step :: Float64
            radius :: Union{Int,Float64}
            noise :: Union{Int,Float64}
            λ₀ :: Union{Int,Float64}
            λ₁ :: Union{Int,Float64}
            λ₂ :: Union{Int,Float64}
            λ₃ :: Union{Int,Float64}
            kill_zone_angle :: Float64 # (true,angle) or (false,compute angle length from γ)
            file_name_path :: String
            iterations :: Int
            random_walk_first_angle::Union{Bool,Float64}
            model::AbstractT6SSMachineModel
        end

        mutable struct MachineStatePosition
            state :: ctm.StateSpace
            position :: CartesianCoordinates
        end
    // # end

    ######################################################
    # StateDependentRandomWalkParameters
    # Get and set functions
    ######################################################
    // # functions
        function get_end_time(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.end_time
        end

        function set_end_time!(sdrw_params::StateDependentRandomWalkParameters,end_time)
            sdrw_params.end_time = Int(end_time)
        end

        function get_time_step(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.time_step
        end


        function set_time_step!(sdrw_params::StateDependentRandomWalkParameters,time_step)
            sdrw_params.time_step = time_step
        end

        function get_radius(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.radius
        end

        function set_radius!(sdrw_params::StateDependentRandomWalkParameters,radius)
            sdrw_params.radius = radius
        end

        function get_noise(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.noise
        end

        function set_noise!(sdrw_params::StateDependentRandomWalkParameters,noise)
            sdrw_params.noise = noise
        end

        function get_λ₀(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.λ₀
        end

        function set_λ₀!(sdrw_params::StateDependentRandomWalkParameters,λ₀)
            sdrw_params.λ₀ = λ₀
        end

        function get_λ₁(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.λ₁
        end

        function set_λ₁!(sdrw_params::StateDependentRandomWalkParameters,λ₁)
            sdrw_params.λ₁ = λ₁
        end

        function get_λ₂(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.λ₂
        end

        function set_λ₂!(sdrw_params::StateDependentRandomWalkParameters,λ₂)
            sdrw_params.λ₂ = λ₂
        end

        function get_λ₃(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.λ₃
        end

        function set_λ₃!(sdrw_params::StateDependentRandomWalkParameters,λ₃)
            sdrw_params.λ₃ = λ₃
        end

        function get_γ(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.γ
        end

        function set_γ!(sdrw_params::StateDependentRandomWalkParameters,γ)
            sdrw_params.γ = γ
        end

        function get_kill_zone_angle(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.kill_zone_angle
        end

        function set_kill_zone_angle!(sdrw_params::StateDependentRandomWalkParameters,kill_zone_angle)
            sdrw_params.kill_zone_angle = kill_zone_angle
        end

        function get_file_name_path(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.file_name_path
        end

        function set_file_name_path!(sdrw_params::StateDependentRandomWalkParameters,file_name_path)
            sdrw_params.file_name_path = file_name_path
        end

        function get_iterations(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.iterations
        end

        function set_iterations!(sdrw_params::StateDependentRandomWalkParameters,iterations)
            sdrw_params.iterations = iterations
        end

        function get_random_walk_first_angle(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.random_walk_first_angle
        end

        function set_random_walk_first_angle!(sdrw_params::StateDependentRandomWalkParameters,random_walk_first_angle)
            sdrw_params.random_walk_first_angle = random_walk_first_angle
        end

        function get_model(sdrw_params::StateDependentRandomWalkParameters)
            sdrw_params.model
        end

        function set_model!(sdrw_params::StateDependentRandomWalkParameters,model::AbstractT6SSMachineModel)
            sdrw_params.model = model
        end
    // # end


    ######################################################
    # Non specific helper functions
    ######################################################
    // # Functions
        rta(t) = round(t,digits = 3)

        function get_ratio_in_each_state(series::Vector{Int64})
            states = unique(series)
            cs = (cumsum([[count(==(i), series[t]) for i in states] for t in eachindex(series)]))
            sum_per_index = [sum(cs[i]) for i in eachindex(series)]
            ratio_per_state = cs ./ sum_per_index
            [getindex.(ratio_per_state, i) for i in states]
        end

        function convert_vector_vector_to_single_vector(series::Vector{Vector{Float64}})
            iterator = eachindex(series[1])
            [getindex.(series, i) for i in iterator]
        end

        function kill_angle(γ::Float64)
            Angle(acos((5-γ^2)/4))
        end

        function angle_length(min_angle,max_angle)
            (max_angle - min_angle)
        end 

        function append_jld2(fname, varname, data)
            jldopen(fname, "a+") do file
            file[varname] = data
            end
        end

        function string_date()
            formatted_today = Dates.format(now(), "yyyymmdd")
            formatted_today
        end


        function file_name_saver(root_path::String,file_details_path::String,extension::String)
            mod_root_path = root_path[end] == '/' ? root_path : string(root_path,"/")
            string(mod_root_path,string_date(),"_",file_details_path,".",lowercase(extension))
        end


    // # end



    ######################################################
    # StateDependentRandomWalkParameters and functions
    ######################################################
    //  # Functions 

        # Waiting time to kill zone functions
        # IC = π ∨ random
        # machines = 1..M
        # λ₁ = ϵ..Λ₁ → π₁ ∈ 0..1
        # n = 1..N (noise)

        function time_series(sdrw_params::StateDependentRandomWalkParameters)
            Δt = get_time_step(sdrw_params)
            T = get_end_time(sdrw_params)
            range(0,T,step = Δt)
        end

        function get_rw_IC_random(sdrw_params::StateDependentRandomWalkParameters)
            rw_angle_bool = get_random_walk_first_angle(sdrw_params)
            boolean_true = rw_angle_bool isa Bool && rw_angle_bool
            boolean_true ? true : false
        end

        function get_num_machines(sdrw_params::StateDependentRandomWalkParameters)
            get_λ₀(sdrw_params)
        end

        function get_state_1_rate(sdrw_params::StateDependentRandomWalkParameters)
            get_λ₁(sdrw_params)
        end

        function get_state_2_rate(sdrw_params::StateDependentRandomWalkParameters)
            get_λ₂(sdrw_params)
        end

        function get_state_3_rate(sdrw_params::StateDependentRandomWalkParameters)
            get_λ₃(sdrw_params)
        end


        function get_gamma(sdrw_params::StateDependentRandomWalkParameters)
            get_γ(sdrw_params)
        end

        function mean_square_displacement_slope(sdrw_params::StateDependentRandomWalkParameters)
            r = get_radius(sdrw_params)
            n = get_noise(sdrw_params)
            (n/r)^2
        end

        function diffusion(sdrw_params::StateDependentRandomWalkParameters)
            msd_slope = mean_square_displacement_slope(sdrw_params)
            (1/2)*msd_slope
        end

        function zero_to_τ(sdrw_params::StateDependentRandomWalkParameters)
            γ = get_gamma(sdrw_params)
            θ(kill_angle(γ))
        end

    // # end

    ####################
    # Time to kill zone
    ####################
    // # functions    
        function get_minimum_angle_length(sdrw_params::StateDependentRandomWalkParameters)
            get_kill_zone_angle(sdrw_params)
        end

        function get_expected_number_cycles(sdrw_params::StateDependentRandomWalkParameters)
            π/get_kill_zone_angle(sdrw_params)
        end



        function get_expected_number_cycles_time_to_state_1(sdrw_params::StateDependentRandomWalkParameters)
            get_expected_number_cycles(sdrw_params)*get_expected_time_to_state_1(sdrw_params)
        end

        function get_angle_length(sdrw_params::StateDependentRandomWalkParameters)
            min_angle = get_minimum_angle_length(sdrw_params)
            max_angle = π
            angle_length(min_angle,max_angle)
        end

        function expected_time_travel_angle_random_IC(sdrw_params::StateDependentRandomWalkParameters)
            L² = get_angle_length(sdrw_params)^2
            D = diffusion(sdrw_params)
            (1/3)*L²/D
        end

        function expected_time_travel_angle_π_IC(sdrw_params::StateDependentRandomWalkParameters)
            L² = get_angle_length(sdrw_params)^2
            D = diffusion(sdrw_params)
            L²/(2*D)
        end
    

        function expected_time_travel(sdrw_params::StateDependentRandomWalkParameters)
            get_rw_IC_random(sdrw_params) ? 
                expected_time_travel_angle_random_IC(sdrw_params) :
                    expected_time_travel_angle_π_IC(sdrw_params)
        end
        
    // # end






    ########################
    # CTMC - State system
    ########################
    // # functions
        Exponential(::State_1,sdrw_params::StateDependentRandomWalkParameters) = Distributions.Exponential(1/get_state_1_rate(sdrw_params))
        Exponential(::State_2,sdrw_params::StateDependentRandomWalkParameters) = Distributions.Exponential(1/get_state_2_rate(sdrw_params))
        Exponential(::State_3,sdrw_params::StateDependentRandomWalkParameters) = Distributions.Exponential(1/get_state_3_rate(sdrw_params))


        function state_markov_simulation_parameters(sdrw_params::StateDependentRandomWalkParameters)

            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            ctm.SimulationParameters(ts,st,tr,π̃)
        end


        function get_probability_per_state(sdrw_params::StateDependentRandomWalkParameters)
            sp = state_markov_simulation_parameters(sdrw_params)
            ctm.get_prob_per_state(sp)
        end

        function get_prbability_state_1(::State_1,sdrw_params::StateDependentRandomWalkParameters)
            get_probability_per_state(sdrw_params)[1]
        end

        function get_prbability_state_2(::State_2,sdrw_params::StateDependentRandomWalkParameters)
            get_probability_per_state(sdrw_params)[1]
        end

        function get_prbability_state_3(::State_3,sdrw_params::StateDependentRandomWalkParameters)
            get_probability_per_state(sdrw_params)[1]
        end
        #=
        function get_expected_time_to_state_1(sdrw_params::StateDependentRandomWalkParameters)
            λ₁ = get_state_1_rate(sdrw_params)
            λ₂ = get_state_2_rate(sdrw_params)
            λ₃ = get_state_3_rate(sdrw_params)
            (1/λ₁ + 1/λ₂ + 1/λ₃)
        end
        =#


        function get_expected_time_to_state_1(sdrw_params::StateDependentRandomWalkParameters)
            λ₁ = get_state_1_rate(sdrw_params)
            λ₂ = get_state_2_rate(sdrw_params)
            λ₃ = get_state_3_rate(sdrw_params)
            rates = [λ₁,λ₂,λ₃]
            π = get_probability_per_state(sdrw_params)
            starting_s1_time = sum(1 ./ rates)
            starting_s2_time = sum(1 ./ rates[2:3])
            starting_s3_time = sum(1 ./ rates[3])
            rates_times = [starting_s1_time,starting_s2_time,starting_s3_time]
            sum(π .* rates_times)
        end

    // # end




    ################################################
    # Random walk and state dependent functions
    ################################################
    // # functions
        #=
        function in_kill_zone(sdrw_params::StateDependentRandomWalkParameters,point::CartesianCoordinates)::Bool
            x,y = coordinates(point)
            kill_angle = Angle(get_kill_zone_angle(sdrw_params))
            rad = Radius(get_radius(sdrw_params))
            p₁ = PolarCoordinates(rad,-kill_angle)
            p₂ = PolarCoordinates(rad,kill_angle)
            c₁ = coordinates(p2c(p₁))
            c₂ = coordinates(p2c(p₂))
            x_kill_zone = (c₁[1] + c₂[1])/2 
            y_kill_zone = (c₁[2],c₂[2])
            x >= x_kill_zone && y_kill_zone[1] <= y <= y_kill_zone[2]
        end
        =#

        function in_kill_zone(sdrw_params::StateDependentRandomWalkParameters,point::CartesianCoordinates)::Bool
            polar_point = c2p(point)
            kill_angle = Angle(get_kill_zone_angle(sdrw_params))
            return abs(kill_angle.θ) > abs(polar_point.θ.θ)
        end
        

        function get_random_angle_or_π(sdrw_params::StateDependentRandomWalkParameters)
            get_rw_IC_random(sdrw_params) ? random_θ() : 1.0*π
        end
        function get_random_walk_init(sdrw_params::StateDependentRandomWalkParameters)
            angle = get_random_angle_or_π(sdrw_params) |> Angle
            radius = Radius(get_radius(sdrw_params))
            p2c(PolarCoordinates(radius,angle))
        end

        function scale_with_state(state::State_1,cn::CircleNoise,s::ScaleNoise)
            CircleNoise(cn.radius,s)
        end
        function scale_with_state(state::Union{State_2,State_3},cn::CircleNoise,s::ScaleNoise)
            CircleNoise(cn.radius,ScaleNoise(0.0))
        end

    // # end




    ################################################
    # State dependent random walk estimates
    ################################################
    // # functions

        function get_state_dependent_coefficient(sdrw_params::StateDependentRandomWalkParameters)
            λ₁ = get_state_1_rate(sdrw_params)
            λ₂ = get_state_2_rate(sdrw_params)
            λ₃ = get_state_3_rate(sdrw_params)
            (1 + λ₁/λ₂ + λ₁/λ₃)
        end



        function get_estimated_waiting_time(::AbstractT6SSMachineModel, sdrw_params::StateDependentRandomWalkParameters) end


        function get_estimated_waiting_time(::RandomWalkOnly, sdrw_params::StateDependentRandomWalkParameters)
            expected_time_travel(sdrw_params)
        end


        function get_estimated_waiting_time(::StateDependentRandomWalk, sdrw_params::StateDependentRandomWalkParameters)
            est_time = get_estimated_waiting_time(RandomWalkOnly(), sdrw_params)
            time_to_be_in_kill_zone_in_right_state = get_expected_number_cycles_time_to_state_1(sdrw_params)
            est_time + time_to_be_in_kill_zone_in_right_state
        end


        function get_estimated_waiting_time(::MultipleMachinesRandomWalkOnly, sdrw_params::StateDependentRandomWalkParameters)
            est_time = get_estimated_waiting_time(RandomWalkOnly(), sdrw_params)
            num_machines = get_num_machines(sdrw_params)
            est_time/num_machines
        end

        function get_estimated_waiting_time(::MultipleMachinesStateDependentRandomWalk, sdrw_params::StateDependentRandomWalkParameters)
            time_to_be_in_kill_zone_in_right_state = get_estimated_waiting_time(StateDependentRandomWalk(), sdrw_params)
            num_machines = get_num_machines(sdrw_params)
            time_to_be_in_kill_zone_in_right_state/num_machines

        end


        function get_estimated_waiting_time(sdrw_params::StateDependentRandomWalkParameters)
            model = get_model(sdrw_params)
            get_estimated_waiting_time(model, sdrw_params)
        end

    // # end



    ################################################
    # Simulations
    ################################################
    // # functions
        
        function get_CTMC(sdrw_params::StateDependentRandomWalkParameters)
            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)
            state = sp.states.s₁
            series = Int[]
            T = ctm.time_series(sp)
            for t in T
                state = ctm.update_single_state(state,tr,ts)
                push!(series,ctm.get_state_int(state))
            end

            series

        end


        function get_simulated_waiting_time(::AbstractT6SSMachineModel, sdrw_params::StateDependentRandomWalkParameters) end

        function get_simulated_waiting_time(::RandomWalkOnly,sdrw_params::StateDependentRandomWalkParameters)
            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)


            # Determine kill zone in XxY
            ikz(step::CartesianCoordinates) = in_kill_zone(sdrw_params,step)
            time_kill_zone = Float64[]
            for n in Base.OneTo(sdrw_params.iterations)
                r = cn.radius
                step = get_random_walk_init(sdrw_params)
                T = time_series(ts_rw)
                for t in T
                    if ikz(step)
                        push!(time_kill_zone,t)
                        break
                    end
                    step = step_along_cartesian(step,cn,ts_rw)

                end
            end
            time_kill_zone
        end

        function get_simulated_waiting_time(::RandomWalkOnlyPerRandomDrawλ₁,sdrw_params::StateDependentRandomWalkParameters)
            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)


            # Determine kill zone in XxY
            ikz(step::CartesianCoordinates) = in_kill_zone(sdrw_params,step)

            time_kill_zone = Float64[]
            num_times_times_walk = Int[]
            for n in Base.OneTo(sdrw_params.iterations)
                r = cn.radius
                state = sp.states.s₁
                step =get_random_walk_init(sdrw_params)
                counter = 1
                temp_time = Float64[]
                for iters in Base.OneTo(1000) # so that num iters till kill zone
                    T = rand(Exponential(State_1(),sdrw_params)) |> t -> range(0,stop = t,step = ts_rw.time_step)
                    for t in T
                        # Check if in kill zone
                        if ikz(step)  
                            push!(time_kill_zone,sum(temp_time) + t)
                            temp_time = Float64[]
                            push!(num_times_times_walk,counter)
                            counter = 1
                            break
                        elseif !ikz(step) && t == T[end]
                            push!(temp_time,t)
                        end
                        step = step_along_cartesian(step,cn,ts_rw)
                    end
                    ikz(step) ? break : counter += 1
                end
            end
            (t_kz = time_kill_zone,n_rws = num_times_times_walk)
        end

        function get_simulated_waiting_time(::StateDependentRandomWalk, sdrw_params::StateDependentRandomWalkParameters)
            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)

            @info "Starting simulation"
            # Determine kill zone in XxY
            ikz(step::CartesianCoordinates) = in_kill_zone(sdrw_params,step)
            time_kill_zone = Float64[]
            for n in Base.OneTo(sdrw_params.iterations)
                r = cn.radius
                state = sp.states.s₁
                step =get_random_walk_init(sdrw_params)
                T = time_series(ts_rw)
                for t in T
                    if ikz(step) && state == sp.states.s₃ 
                        push!(time_kill_zone,t)
                        break
                    elseif !ikz(step) && t == T[end]
                        error("Not in kill zone and sim ended")
                    end
                    state = ctm.update_single_state(state,tr,ts)
                    cn = scale_with_state(state,cn,s)
                    step = step_along_cartesian(step,cn,ts_rw)
                end
            end
            time_kill_zone
        end

        function get_simulated_waiting_time(::MultipleMachinesRandomWalkOnly, sdrw_params::StateDependentRandomWalkParameters)
        end

        function get_state_dependent_noise(sdrw_params::StateDependentRandomWalkParameters,state)
            r = Radius(get_radius(sdrw_params))
            s = ScaleNoise(get_noise(sdrw_params))
            cn = CircleNoise(r,s)
            scale_with_state(state,cn,s)
        end

        function get_random_walk_module_TimeSeries_struct(sdrw_params::StateDependentRandomWalkParameters)
            rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
        end
        
        function update_state(sp,state)
            ctm.update_single_state(state,sp.transition_rates,sp.time_series)
        end



        function update_process(
            ::Union{StateDependentRandomWalk,MultipleMachinesStateDependentRandomWalk},
            sdrw_params::StateDependentRandomWalkParameters,process_param)
            state,step = process_param
            state = update_state(sp,state)
            rw_based_ts = get_random_walk_module_TimeSeries_struct(sdrw_params)
            scaled_noise = get_state_dependent_noise(sdrw_params,state)
            step = step_along_cartesian(step,scaled_noise,rw_based_ts)
            (state,step)
        end

        function update_process(
            ::Union{RandomWalkOnly,MultipleMachinesRandomWalkOnly,RandomWalkOnlyPerRandomDrawλ₁},
            sdrw_params::StateDependentRandomWalkParameters,process_param)
            step = process_param
            state = State_1()
            rw_based_ts = get_random_walk_module_TimeSeries_struct(sdrw_params)
            scaled_noise = get_state_dependent_noise(sdrw_params,state)
            step = step_along_cartesian(step,scaled_noise,rw_based_ts)
            step
        end





        function get_simulated_waiting_time(::MultipleMachinesStateDependentRandomWalk, sdrw_params::StateDependentRandomWalkParameters)
            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            weights = get_probability_per_state(sdrw_params)
            state = sample([ctm.State_1(),ctm.State_2(),ctm.State_3()],Weights(weights))
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)

            N_mach = sdrw_params.λ₀#get_num_machines(sdrw_params.λ₀)

            # Determine kill zone in XxY
            ikz(step::CartesianCoordinates) = in_kill_zone(sdrw_params,step)
            time_kill_zone = Float64[]
            for n in Base.OneTo(sdrw_params.iterations)
                machines_min_kill_time = Float64[]
                for m in Base.OneTo(N_mach)
                    r = cn.radius
                    state = sample([ctm.State_1(),ctm.State_2(),ctm.State_3()],Weights(weights))
                    step =get_random_walk_init(sdrw_params)
                    T = time_series(ts_rw)
                    for t in T
                        if ikz(step) && state == sp.states.s₃
                            push!(machines_min_kill_time,t)
                            break
                        elseif !ikz(step) && t == T[end]
                            error("Not in kill zone and sim ended")
                        end
                        state = ctm.update_single_state(state,tr,ts)
                        cn = scale_with_state(state,cn,s)
                        step = step_along_cartesian(step,cn,ts_rw)
                        # Check if in kill zone
                    end
                end
                push!(time_kill_zone,minimum(machines_min_kill_time))
            end
            time_kill_zone
        end

        function get_simulated_waiting_time(::StateIndependentRandomWalk, sdrw_params::StateDependentRandomWalkParameters)
            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)

            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)


            # Determine kill zone in XxY
            ikz(step::CartesianCoordinates) = in_kill_zone(sdrw_params,step)

            time_kill_zone = Float64[]
            
            for n in Base.OneTo(sdrw_params.iterations)
                r = cn.radius
                state = sp.states.s₁
                step =get_random_walk_init(sdrw_params)
                T = time_series(ts_rw)

                for t in T
                    # Check if in kill zone
                    if ctm.get_state_int(state) == 3 && ikz(step)
                        push!(time_kill_zone,t)
                        break
                    end
                    state = ctm.update_single_state(state,tr,ts)
                    step = step_along_cartesian(step,cn,ts_rw)
                end
            end
            time_kill_zone
        end

        function get_simulated_waiting_time(sdrw_params::StateDependentRandomWalkParameters)
            model = get_model(sdrw_params)
            get_simulated_waiting_time(model, sdrw_params)
        end

    

    // # end


    ################################################
    # Animation structs and functions
    ################################################
    // # functions

        function run_simulation_save_animation(sdrw_params::StateDependentRandomWalkParameters)
            # Animate state dependent random walk on a circle 
            # Colors
            colormap = Dict(1 => RGB(0, 1, 0),  # Green
                            2 => RGB(1, 0.5, 0),  # Orange
                            3 => RGB(1, 0, 0))  # Red
            # Random walk parameters

            ts_rw = rw.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            r = Radius(sdrw_params.radius)
            s = ScaleNoise(sdrw_params.noise)
            cn = CircleNoise(r,s)
            state = State_1()
            cn = scale_with_state(state,cn,s)





            # Markov process parameters
            ts = ctm.TimeSeries(sdrw_params.end_time,sdrw_params.time_step)
            st = ctm.States(ctm.State_1(),ctm.State_2(),ctm.State_3())
            tr = ctm.Rates(sdrw_params.λ₁,sdrw_params.λ₂,sdrw_params.λ₃)
            π̃ = ctm.Π([1/3,1/3,1/3])
            sp = ctm.SimulationParameters(ts,st,tr,π̃)



            r = cn.radius
            state = sp.states.s₁
            step = random_cartesian_circle_point(r)
            state_step = MachineStatePosition[]
            push!(state_step, MachineStatePosition(state, step))
            T = time_series(ts_rw)
            for t in T
                state = ctm.update_single_state(state,tr,ts)
                cn = scale_with_state(state,cn,s)
                step = step_along_cartesian(step,cn,ts_rw)
                push!(state_step,MachineStatePosition(state, step))
            end


            # Set up figure
            points = gl.Observable(cart2Point2f(state_step[1].position))
            machine_state = gl.Observable([colormap[ctm.get_state_int(state_step[1].state)]])

            r = Radius(sdrw_params.radius)
            ϕ₁=Angle(get_kill_zone_angle(sdrw_params))
            p₁ = PolarCoordinates(r,ϕ₁)
            p₂ = PolarCoordinates(r,-ϕ₁)
            (x₁,y₁) = coordinates(p2c(p₁))
            (x₂,y₂) = coordinates(p2c(p₂))
            x₁_range = range(0,x₁,length = 100)
            x₂_range = range(0,x₂,length = 100)
            y₁_range = range(0,y₁,length = 100)
            y₂_range = range(0,y₂,length = 100)

            # Define the angles for the arc
            angles = range(θ(-ϕ₁), stop=θ(ϕ₁), length=100)

            # Create the points for the polygon
            polygon_points = [Point2f(0,0), Point2f(x₂,y₂)]
            append!(polygon_points, [Point2f(radius(r)*cos(θ), radius(r)*sin(θ)) for θ in angles])
            push!(polygon_points, Point2f(x₁,y₁))

            # Set up circle
            ϕ = gl.LinRange(0, 2π, 100)
            ps = PolarCoordinates.(Ref(r),Angle.(ϕ))
            cs = coordinates.(p2c.(ps))

            # Set up figure
            fig = gl.Figure(size = (800, 800),fontsize = 25)
            ax = gl.Axis(fig[1, 1],aspect = 1,xlabel="x",ylabel="y")
            gl.lines!(ax,cs,label = "Cell",color = :orange)
            gl.lines!(ax,x₁_range,y₁_range,color = :red)
            gl.lines!(ax,x₂_range,y₂_range,color = :red)
            gl.poly!(ax, polygon_points, color = (:red, 0.4)) 
            gl.arc!(ax,Point2f(0,0),radius(r),θ(-ϕ₁),θ(ϕ₁),color = :red,label = "Kill zone")
            gl.scatter!(ax,points,markersize = 30,color = machine_state)
            gl.axislegend()


            frame_rate = 60
            gl.record(fig, sdrw_params.file_name_path, state_step;
                    framerate = frame_rate) do s
            points[] = cart2Point2f(s.position)
            machine_state[] = [colormap[ctm.get_state_int(s.state)]]
            end


            isfile(sdrw_params.file_name_path) ? println("Success! File exists") : println("Failure! File does not exist")
        end
        
        function save_single_random_walk_on_circle_video(ts::rw.TimeSeries,cn::CircleNoise,file_name_path,frame_rate)
            r = cn.radius
        
        
            T = time_series(ts)
            step = random_cartesian_circle_point(r)
            steps = []
            for t in T
                step = step_along_cartesian(step,cn,ts)
                push!(steps,step)
            end


            # Set up figure
            points = gl.Observable(cart2Point2f(steps[1]))
            machine_state



            ϕ = gl.LinRange(0, 2π, 100)
            ps = PolarCoordinates.(Ref(r),Angle.(ϕ))
            cs = coordinates.(p2c.(ps))
            fig = gl.Figure(size = (800, 800),fontsize = 25)
            ax = gl.Axis(fig[1, 1],aspect = 1,xlabel="x",ylabel="y")
            gl.lines!(ax,cs,label = "Cell",color = :orange)
            gl.scatter!(ax,points,markersize = 30,color = machine_state)
            gl.axislegend()
            fig





            # Save animation to file 
            gl.record(fig, file_name_path, steps;
                    framerate = frame_rate) do s
                points[] = cart2Point2f(s)
            end

            isfile(file_name_path) ? println("Success! File exists") : println("Failure! File does not exist")
        end
    // # end



end