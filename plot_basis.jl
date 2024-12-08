using Plots

# STO와 GTO 기저 함수 정의
function sto(r, n, etha)
    return r^(n-1) * exp(-etha*r)
end

function gto(r, n, alpha)
    return r^(n-1) * exp(-alpha*r^2)
end

# 비교 플롯 생성
function compare_basis()
    r = 0:0.1:10
    
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
    
    plot(p1, p2, layout=(2,1), size=(800,600))
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
    r = 0:0.1:10
    
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
