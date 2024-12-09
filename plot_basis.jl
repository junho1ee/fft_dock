using Plots

# STO와 GTO 기저 함수 정의
function sto(r, n, etha)
    return abs(r)^(n-1) * exp(-etha*abs(r))
end

function gto(r, n, alpha)
    return abs(r)^(n-1) * exp(-alpha*r^2)
end

# 비교 플롯 생성
function compare_basis()
    r = -5:0.1:5
    
    # 1s orbital (n=1)
    sto_1s = sto.(r, 1, 1.0)
    gto_1s = gto.(r, 1, 0.5)  # alpha 값은 대략적으로 조정됨
    
    # 2s orbital (n=2)
    sto_2s = sto.(r, 2, 1.0)
    gto_2s = gto.(r, 2, 0.5)
    
    p1 = plot(r, [sto_1s gto_1s], 
        label=["1s STO" "1s GTO"],
        title="1s orbital comparison",
        xlabel="distance (Bohr)",
        ylabel="amplitude")
    
    p2 = plot(r, [sto_2s gto_2s], 
        label=["2s STO" "2s GTO"],
        title="2s orbital comparison",
        xlabel="distance (Bohr)",
        ylabel="amplitude")
    
    plot(p1, p2, layout=(2,1), size=(800,800))
end

# STO-nG 기저함수 정의
function sto_ng(r, n, ng_params)
    result = 0.0
    for (coef, alpha) in ng_params
        result += coef * gto(r, n, alpha)
    end
    return result
end

# STO-3G와 STO-6G 매개변수 (1s 오비탈)
sto3g_params = [
    (0.444635, 0.109818),
    (0.535328, 0.405771),
    (0.154329, 2.227660)
]

sto6g_params = [
    (0.195233, 0.050331),
    (0.405345, 0.151399),
    (0.349938, 0.384597),
    (0.109391, 1.233787),
    (0.015432, 4.683033),
    (0.001523, 20.638684)
]

function compare_basis_extended()
    r = -5:0.1:5
    
    # 1s 오비탈 비교
    sto_1s = sto.(r, 1, 1.0)
    gto_1s = gto.(r, 1, 0.5)
    sto3g_1s = [sto_ng(ri, 1, sto3g_params) for ri in r]
    sto6g_1s = [sto_ng(ri, 1, sto6g_params) for ri in r]
    
    p = plot(r, [sto_1s gto_1s sto3g_1s sto6g_1s], 
        label=["STO" "single GTO" "STO-3G" "STO-6G"],
        title="1s orbital comparison",
        xlabel="distance (Bohr)",
        ylabel="amplitude",
        legend=:topright)
    
    # 차이 플롯
    p_diff = plot(r, [abs.(sto3g_1s .- sto_1s) abs.(sto6g_1s .- sto_1s)], 
        label=["STO-3G error" "STO-6G error"],
        title="difference from STO",
        xlabel="distance (Bohr)",
        ylabel="absolute error",
        yscale=:log10)
    
    plot(p, layout=(1,1), size=(800,400))
end

# 확장된 비교 플롯 생성
compare_basis_extended() 

# 플롯 생성 및 표시
# compare_basis()

# 구면 조화 함수 (l=2, m=0 for 3d_z2)
function Y_20(theta, phi)
    return 0.25 * sqrt(5/π) * (3*cos(theta)^2 - 1)
end

# 3d 오비탈 함수 (3d_z2)
function orbital_3d(r, theta, phi, type="sto")
    # 반지름 부분
    if type == "sto"
        R = r^2 * exp(-r)  # n=3 for 3d orbital
    else  # gto
        R = r^2 * exp(-0.3*r^2)  # alpha 값은 적절히 조정
    end
    
    # 각도 부분 (3d_z2 오비탈)
    Y = Y_20(theta, phi)
    
    return R * Y
end

function plot_3d_orbital()
    # 구면 좌표계 그리드 생성
    r = range(0, 5, length=50)
    theta = range(0, π, length=50)
    phi = range(0, 2π, length=50)
    
    # 메쉬그리드 생성
    R = [r_i for r_i in r, _ in theta]
    Θ = [theta_i for _ in r, theta_i in theta]
    
    # 카테시안 좌표로 변환
    X = R .* sin.(Θ)
    Y = zeros(size(X))  # phi=0 평면에서의 단면
    Z = R .* cos.(Θ)
    
    # 오비탈 값 계산
    V_sto = [orbital_3d(R[i,j], Θ[i,j], 0.0, "sto") for i in 1:size(R,1), j in 1:size(R,2)]
    V_gto = [orbital_3d(R[i,j], Θ[i,j], 0.0, "gto") for i in 1:size(R,1), j in 1:size(R,2)]
    
    # 플롯 생성
    p1 = surface(X, Z, V_sto,
        title="3d_z² STO orbital",
        xlabel="x",
        ylabel="z",
        zlabel="ψ",
        camera=(45, 45))
    
    p2 = surface(X, Z, V_gto,
        title="3d_z² GTO orbital",
        xlabel="x",
        ylabel="z",
        zlabel="ψ",
        camera=(45, 45))
    
    plot(p1, p2, layout=(1,2), size=(1200,500))
end

# 3D 오비탈 플롯 생성
plot_3d_orbital()
