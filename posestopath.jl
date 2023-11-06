using LinearAlgebra
using StaticArrays, Rotations
using Plots

mutable struct SE2
    R::Float64
    p::Vector{Float64}
    SE2(theta::Float64, p) = new(theta, p)
    SE2(R::Matrix{Float64}, p) = new(rotation_angle(R), p)
    SE2(R::RotMatrix2{Float64}, p) = new(rotation_angle(R), p)
end

mutable struct se2
    omega::Float64
    v::Vector{Float64}
end

function wedge(a::Float64)
    return [
        0.0 -a;
        a    0.0
    ]    
end

function logvee(pose::SE2)
    omega = pose.R
    if abs(omega) < 1e-9
        v = pose.p
    else
        v = wedge(-omega) * inv(I - RotMatrix{2,Float64}(omega)) * pose.p
    end
    return se2(omega, v)
end

function expwedge(twist::se2, T)
    R = twist.omega * T
    if abs(R) < 1e-9
        p = twist.v * T
    else
        p = (I - RotMatrix{2,Float64}(R)) * wedge(1.0) * twist.v / twist.omega
    end
    return SE2(R, p)
end

function Pose2Matrix(pose::SE2)
    return [RotMatrix{2,Float64}(pose.R) pose.p; zeros(1, 2) 1]
end
function Matrix2Pose(M) 
    return SE2(rotation_angle(RotMatrix{2,Float64}(M[1:2, 1:2])), M[1:2, 3])
end

function SolveOptimalPath(pose_dst::SE2)
    twist_dst = logvee(pose_dst)

    aligned_twist = se2(0.0, [0.0; 0.0])
    aligned_twist.omega = twist_dst.omega
    aligned_twist.v[1] = sign(twist_dst.v[1]) * norm(twist_dst.v)

    required_lateral_displacement = atan(twist_dst.v[2], abs(twist_dst.v[1])) * norm(twist_dst.v)

    Y_f = required_lateral_displacement
    X_f = [aligned_twist.omega, aligned_twist.v[1]]
        
    K = 4 * abs(Y_f) / norm(X_f)^2

    function FuncForBisection(theta)
        return (theta - sin(theta)) / (1 - cos(theta))
    end

    left, right = 1e-9, 2pi
    mid = 0.0
    for i = 1:25
        mid = (left + right) / 2.0
        if (FuncForBisection(mid) < K)
            left = mid
        else
            right = mid
        end
    end

    theta_puseudo = sign(Y_f) * mid
    twist_puseudo = logvee(SE2(theta_puseudo, X_f))
    return twist_puseudo
end

scale = 5.0 #relative weight of translation againt the rotation. recommend to use minimum arc lenght of the vehicle as the weight
x_f = 45 * scale
y_f = 10 * scale

pose_src_fin = SE2(0.0, [0.0, 0.0])
pose_dst_fin = SE2(-0.7, [x_f, y_f])

pose_diff = Matrix2Pose(Pose2Matrix(pose_dst_fin) * inv(Pose2Matrix(pose_src_fin)))
pose_diff_modified = Matrix2Pose(Pose2Matrix(pose_dst_fin) * inv(Pose2Matrix(pose_src_fin)))
error_rate_of_result_point = 1e9

trajectories = []

for num_trial = 1:20
    @time twist_puseudo = SolveOptimalPath(pose_diff_modified)
    pose_temp = Pose2Matrix(SE2(0.0, [0.0, 0.0]))
    dt = 0.01
    x_traj, y_traj = [], []
    @time for time_puseudo = range(0, 1.0, step = dt)
        input = RotMatrix{2,Float64}(twist_puseudo.omega * time_puseudo) * twist_puseudo.v
        diff = Pose2Matrix(expwedge(se2(input[1], [input[2], 0.0]), dt))
        pose_temp = pose_temp * diff
    
        push!(x_traj, pose_temp[1, 3])
        push!(y_traj, pose_temp[2, 3])      
    end
    push!(trajectories, (x_traj, y_traj))

    pose_f = Matrix2Pose(pose_temp)
    error_of_point = pose_f.p - pose_diff.p
    error_rate_of_result_point = norm(error_of_point) / norm(pose_diff.p)
 
    if error_rate_of_result_point > 1.0
        # println("Failue Resulting pose is $(pose_f.R), $(pose_f.p)")
        break
    elseif error_rate_of_result_point < 0.01
        # println("Complete to solve. Resulting pose is $(pose_f.R), $(pose_f.p)")
        plot
        break
    else
        # println("Success to solve. Resulting pose of trial $num_trial is $(pose_f.R), $(pose_f.p)")
        pose_diff_modified.p -= error_of_point .* 0.9/(error_rate_of_result_point+1.0)
    end
end

plt = plot()
for traj in trajectories
    plot!(plt, traj ./ scale, aspect_ratio = :equal, lw=2, ls=:dot, size=(800,600))
end
plot(plt)


#path shifterの代わり
#A*のノルム
#数値最適化の基底関数
#数値最適化の初期解


