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

# 플롯 생성 및 표시
compare_basis()
