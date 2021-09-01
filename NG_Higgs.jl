using SpecialFunctions
using QuadGK
using Base.Threads
# using Interpolations
# using DataInterpolations
using DelimitedFiles
using Plots
# plotlyjs()
using PyCall
const itp = pyimport("scipy.interpolate")
const itg = pyimport("scipy.integrate")
const kpivot, ϕpivot=0.05, 1.40
const PcsBack,PcsPivot,PcsIni,PcsPert=0.00025,1e-4,0.0001,0.00025
const RatioMin=1000
const κ=6e-2
const ks,kl=1e12,5e13
const d,c,xc,v0=1.05e10,2.04e-10,1.344,1.24e-9
const ifprint=0
function Rk4(y::Array,n,k,h)
    h2=h/2.
    h6=h/6.
    dy1 = Derivs(y,k,n)
    yt  = y .+ h2*dy1
    dy2 = Derivs(yt,k,n)
    yt  = y .+ h2*dy2
    dy3 = Derivs(yt,k,n)
    yt  = y .+  h*dy3
    dy3 = dy2 .+ dy3
    dy2 = Derivs(yt,k,n)
    return y .+ h6*(dy1 .+ 2dy3 .+ dy2)
end
function Potential(ϕ)
    return [v0*ϕ^(4.0)/4, v0*ϕ^(3.0), 3*v0*ϕ^(2)]
end
function Gfun(x)
    if x>xc
        return [-d/((x-xc)/c+1.)/2-x^(22)/2, d/c/(1+(x-xc)/c)^2/2-11*x^(21), -d/c^2/(1+(x-xc)/c)^3-231*x^(20)]
    else
        return [-d/(-(x-xc)/c+1.)/2-x^(22)/2, -d/c/(1-(x-xc)/c)^2/2-11*x^(21), -d/c^2/(1-(x-xc)/c)^3-231*x^(20)]
    end
end
function Hubble(dϕ,g,V)
    return sqrt((dϕ^2/2+V-g*dϕ^2)/3)
end

function Derivs(y::Array,k,n)
    V,dV,ddV=Potential(y[2])
    g,dg=Gfun(y[2])#dg=gp
    dϕ=y[3]/y[1]
    H=Hubble(dϕ,g,V)#sqrt((1/2*dϕ^2+V-dϕ^2*g)/3)
    ddϕ=-3H*dϕ-(dV-dϕ^2*dg)/(1-2g)
    ϵH=dϕ^2*(1-2g)/H^2/2
    ηH=3+(dV-dϕ^2*dg)/(1-2g)/H/dϕ
    η=2*(ϵH-ηH-dg*dϕ/H/(1-2g))
    zpz=y[1]*H*(1+η/2)
    if n==3
        return [y[1]^2*H,y[3],y[1]^2*(ddϕ+dϕ*H)]
    else
        return [y[1]^2*H,y[3],y[1]^2*(ddϕ+dϕ*H),y[6],y[7],-2zpz*y[6]-k^2*y[4],-2zpz*y[7]-k^2*y[5]]
    end
end
function EvolveEqns(k,yi::Array,accuracy)
    y  = yi
    dy = Derivs(y,k,3)
    aH = dy[1]/y[1]
    dτ = accuracy*min(1/aH,abs(y[3]/dy[3]))
    y  = Rk4(y,3,k,dτ)
    return [y,dτ]
end
function EvolveBackground(yi::Array,ϕs,accuracy)
    k=0
    y  = yi
    dy = Derivs(y,k,3)
    aH = dy[1]/y[1]
    dτ = accuracy*min(1/aH,abs(y[3]/dy[3]))
    V,dV = Potential(y[2])
    if dV<=0
        while y[2]<=(ϕs-y[3]*dτ)
            y,dτ=EvolveEqns(k,y,accuracy)
        end
    else
        while y[2]>(ϕs-y[3]*dτ)
            y,dτ=EvolveEqns(k,y,accuracy)
        end
    end
    dy = Derivs(y,k,3)
    dτ = (ϕs-y[2])/dy[2]
    return y .+ dy*dτ
end
function CheckSR(ϕ,dϕ)
    V,dV,ddV=Potential(ϕ)
    g,dg=Gfun(ϕ)
    H=Hubble(dϕ,g,V)#sqrt((dϕ^2/2+V-dϕ^2*g)/3)
    ϵH=(dϕ^2/2*(1-2g)/H^2)
    ηH=3+(dV-dg*dϕ^2)/(1-2g)/H/dϕ
    η=2*(ϵH-ηH-dg*dϕ/H/(1-2g))
    dη=2*(2H*ϵH^2-3H*ϵH*ηH-H*ηH^2+3H*ϵH+3H*ηH-3*(ηH-2+ϵH)*dg*dϕ/(1-2g)-ddV/(1-2g)/H-2*(dg*dϕ)^2/(1-2g)^2/H)/H
    return [ϵH,η,dη,ηH]
end
function Attractor(ϕ,accuracy,precision)
    maxtime=1000
    V,dV=Potential(ϕ)
    g,=Gfun(ϕ)
    H=sqrt(V/3)
    dϕ= -dV/3H/(1-2g)
    counter=0
    ϕ0=ϕ
    dϕ0=dϕ
    dϕn=dϕ0/1e6
    y=Float64[]
    while abs(dϕ0/dϕn-1)>precision
        counter+=1
        if counter>maxtime
            throw(ErrorException("Couldn't find an attractor after $maxtime times of iterations"))
        end
        dϕn=dϕ0
        ϕ0+= -3dϕ/H
        V1,dV1=Potential(ϕ0)
        g1,=Gfun(ϕ0)
        H=sqrt(V1/3)
        dϕ= -dV/3H/(1-2g1)
        y=[1,ϕ0,dϕ]
        y=EvolveBackground(y,ϕ,accuracy)
        dϕ0=y[3]/y[1]
    end
    dϕ=y[3]/y[1]
    H=Hubble(dϕ,g,V)#sqrt((dϕ^2/2+V-dϕ^2*g)/3)
    ϵH,η,dη,ηH=CheckSR(ϕ,dϕ)
    # display(dϕ)
    ifprint==1 ? display("at ϕ0=$ϕ,dϕ0=$(dϕ/H),ϵH=$ϵH,η=$η,ηH=$ηH") : nothing
    return [dϕ,H]
end
function ReachInf(yi::Array,accuracy,ai)
    k=0
    ϕi=yi[2]
    #display("The initial scale factor a=$ai, the initial ϕi=$ϕi,N=$(log(yi[1]/ai))")
    dy =Derivs(yi,k,3)
    V,=Potential(yi[2])
    g,=Gfun(yi[2])
    dϕ=yi[3]/yi[1]
    H=Hubble(dϕ,g,V)#sqrt((dϕ^2/2+V-g*dϕ^2)/3)
    while dϕ^2*(1-2g)/(2H^2)<1
        dy  = Derivs(yi,k,3)
        aH  = dy[1]/yi[1]
        dτ  = accuracy*min(1/aH,abs(yi[3]/dy[3]))
        yi  = Rk4(yi,3,k,dτ)
        V,  = Potential(yi[2])
        g,  = Gfun(yi[2])
        dϕ  = yi[3]/yi[1]
        H   = Hubble(dϕ,g,V)#sqrt((dϕ^2/2+V-g*dϕ)/3)
    end
    dy  = Derivs(yi,k,3)
    aH  = dy[1]/yi[1]
    ϕe  = yi[2]
    ae  = yi[1]
    ke  = dy[1]/yi[1]
    Ne  = log(yi[1]/ai)
    ϵ   = dϕ^2*(1-2g)/2H^2
    return [ϕe,ke,ae,Ne,ϵ]
end
function ReachAH(yi::Array,aHs,accuracy)
    k   = 0
    y   = yi
    dy  = Derivs(y,k,3)
    aH  = dy[1]/y[1]
    dτ  = accuracy*min(1/aH,abs(y[3]/dy[3]))
    while aH<=aHs
        y   = Rk4(y,3,k,dτ)
        dy  = Derivs(y,k,3)
        aH  = dy[1]/y[1]
        dτ  = accuracy*min(1/aH,abs(y[3]/dy[3]))
    end
    return y
end
function ReachAHBack(y::Array,aHs,accuracy)
    k   = 0
    dy .= Derivs(y,k,3)
    aH  = dy[1]/y[1]
    dτ  = -accuracy*min(1/aH,abs(y[3]/dy[3]))
    while aH>=aHs
        y  .= Rk4(y,3,dτ)
        dy .= Derivs(y,k,3)
        aH  = dy[1]/y[1]
        dτ  = -accuracy*min(1/aH,abs(y[3]/dy[3]))
    end
    return y
end
function EvolveK(yi::Array,k,kmax,accuracy,ae,ai)
    NOs=Float64[]
    NIs=Float64[]
    ReζOs,ImζOs=Float64[],Float64[]
    RedζOs,ImdζOs=Float64[],Float64[]
    ReζIs,ImζIs=Float64[],Float64[]
    RedζIs,ImdζIs=Float64[],Float64[]
    τ=Float64[]
    ϵHOs,ϵHIs = Float64[],Float64[]
    aOs,aIs   = Float64[],Float64[]
    HOs,HIs   = Float64[],Float64[]
    ηOs,ηIs   = Float64[],Float64[]
    dηOs,dηIs = Float64[],Float64[]
    τ0 = 0
    y    = ReachAH(yi,k/RatioMin,PcsBack)
    dy   = Derivs(y,k,3)
    dV   = Potential(y[2])[2]
    H    = dy[1]/y[1]^2
    aH   = H*y[1]
    dϕ   = dy[2]/y[1]
    para = CheckSR(y[2],dϕ)
    ϵH,η = para
    g,   = Gfun(y[2])
    z    = y[1]*dϕ/H*sqrt(1-2g)
    zpz  = y[1]*H*(1+η/2)
    yini = [1/sqrt(2k)/z,0,-zpz/sqrt(2k)/z,-k/sqrt(2k)/z]
    y    = vcat(y,yini)
    dτ   = accuracy*2π/k
    (kmax<=5e17) ? (ratio=1e-5) : (ratio=1e-3)
    while y[1]<ae && kmax/aH>=ratio
        y  = Rk4(y,11,k,dτ)
        dy = Derivs(y,k,11)
        dϕ = dy[2]/y[1]
        para = CheckSR(y[2],dϕ)
        H  = dy[1]/y[1]^2
        aH = H*y[1]
        if (kmax/aH)<=RatioMin
            if kmax/aH>100
                push!(ReζIs,y[4])
                push!(ImζIs,y[5])
                push!(RedζIs,y[6])
                push!(ImdζIs,y[7])
                push!(aIs,y[1])
                push!(HIs,H)
                push!(ϵHIs,para[1])
                push!(ηIs,para[2])
                push!(dηIs,para[3]*aH)
                push!(NIs,log(y[1]/ai))
                τ0 += dτ
                push!(τ,τ0)
            else
                push!(ReζOs,y[4])
                push!(ImζOs,y[5])
                push!(RedζOs,y[6]/aH)#dζ=dζ/dN
                push!(ImdζOs,y[7]/aH)
                push!(NOs,log(y[1]/ai))
                push!(aOs,y[1])
                push!(HOs,H)
                push!(ϵHOs,para[1])
                push!(ηOs,para[2])
                push!(dηOs,para[3])
            end
        end
        dτ = accuracy*2π/max(sqrt(abs(dy[6]/y[4])),k,aH,abs(dy[3]/y[3]))
    end
    NOs    = vcat(last(NIs),NOs)
    ReζOs  = vcat(last(ReζIs),ReζOs)
    ImζOs  = vcat(last(ImζIs),ImζOs)
    ImdζOs = vcat(last(ImdζIs)/last(aIs)/last(HIs),ImdζOs)
    RedζOs = vcat(last(RedζIs)/last(aIs)/last(HIs),RedζOs)
    aOs    = vcat(last(aIs),aOs)
    HOs    = vcat(last(HIs),HOs)
    ϵHOs   = vcat(last(ϵHIs),ϵHOs)
    ηOs    = vcat(last(ηIs),ηOs)
    dηOs   = vcat(last(dηIs)/last(aIs)/last(HIs),dηOs)
    τ      = τ .-τ0
    return [NOs,[ReζIs,ImζIs,RedζIs,ImdζIs],[ReζOs,ImζOs,RedζOs,ImdζOs],[aIs,HIs,ϵHIs,ηIs,dηIs],[aOs,HOs,ϵHOs,ηOs,dηOs],τ]
end
function Kernelτ(a ,H ,ϵH ,dη,k1,k2,k3,Reζτ1 ,Reζτ2 ,Reζτ3 ,Imζτ1 ,Imζτ2 ,Imζτ3 ,Redζτ1 ,Redζτ2 ,Redζτ3 ,Imdζτ1 ,Imdζτ2 ,Imdζτ3 )
    perk1  = complex(Reζτ1,Imζτ1)
    perk2  = complex(Reζτ2,Imζτ2)
    perk3  = complex(Reζτ3,Imζτ3)
    dperk1 = complex(Redζτ1,Imdζτ1)
    dperk2 = complex(Redζτ2,Imdζτ2)
    dperk3 = complex(Redζτ3,Imdζτ3)
    BG     = a^2*ϵH^2
    BG4    = a^2*ϵH*dη
    BG5    = a^2*ϵH^3
    BG6    = BG5
    per1   = conj(perk1)*conj(dperk2)*conj(dperk3)+conj(dperk1)*conj(perk2)*conj(dperk3)+conj(dperk1)*conj(dperk2)*conj(perk3)
    per2   = (-k1^2-k2^2-k3^2)*conj(perk1)*conj(perk2)*conj(perk3)
    per3   = (  ( (k2^2-k1^2-k3^2)/k1^2 + (k1^2-k2^2-k3^2)/k2^2 )*conj(dperk1)*conj(dperk2)*conj(perk3)
                +((k2^2-k3^2-k1^2)/k3^2 + (k3^2-k2^2-k1^2)/k2^2)*conj(dperk3)*conj(dperk2)*conj(perk1)
                +((k3^2-k2^2-k1^2)/k1^2 + (k1^2-k3^2-k2^2)/k3^2)*conj(dperk1)*conj(dperk3)*conj(perk2)  )
    per4   = conj(perk1)*conj(perk2)*conj(dperk3) + conj(dperk1)*conj(perk2)*conj(perk3) + conj(perk1)*conj(dperk2)*conj(perk3)
    per5   = (((k2^2-k1^2-k3^2)/k1^2+(k1^2-k2^2-k3^2)/k2^2 )*conj(dperk1)*conj(dperk2)*conj(perk3)
                +((k2^2-k3^2-k1^2)/k3^2+(k3^2-k2^2-k1^2)/k2^2 )*conj(dperk3)*conj(dperk2)*conj(perk1)
                +((k3^2-k1^2-k2^2)/k1^2+(k1^2-k3^2-k2^2)/k3^2 )*conj(dperk1)*conj(dperk3)*conj(perk2))
    per6   = (k3^2*(k3^2-k1^2-k2^2)/k1^2/k2^2*conj(dperk1)*conj(dperk2)*conj(perk3)
                +k1^2*(k1^2-k3^2-k2^2)/k3^2/k2^2*conj(dperk3)*conj(dperk2)*conj(perk1)
                +k2^2*(k2^2-k1^2-k3^2)/k1^2/k3^2*conj(dperk1)*conj(dperk3)*conj(perk2))
    # temp1  = -4BG*per1
    # temp2  = 2BG*per2
    # temp3  = 2BG*per3
    # temp4  = -2BG4*per4
    # temp5  = -BG5*per5/2
    # temp6  = BG6*per6/2
    return sum([-4BG*per1,2BG*per2,2BG*per3,-2BG4*per4,-BG5*per5/2,-BG6*per6/2])
end
function KernelN(a ,H ,ϵH ,dη,k1,k2,k3,ReζN1 ,ReζN2 ,ReζN3 ,ImζN1 ,ImζN2 ,ImζN3 ,RedζN1 ,RedζN2 ,RedζN3 ,ImdζN1 ,ImdζN2 ,ImdζN3 )
    perk1  = complex(ReζN1,ImζN1)
    perk2  = complex(ReζN2,ImζN2)
    perk3  = complex(ReζN3,ImζN3)
    dperk1 = complex(RedζN1,ImdζN1)
    dperk2 = complex(RedζN2,ImdζN2)
    dperk3 = complex(RedζN3,ImdζN3)
    BG     = a^3*ϵH^2*H
    BG2    = a*ϵH^2/H
    BG3    = BG
    BG4    = a^3*ϵH*dη*H
    BG5    = a^3*ϵH^3*H
    BG6    = BG5
    per1   = conj(perk1)*conj(dperk2)*conj(dperk3) + conj(dperk1)*conj(dperk2)*conj(perk3) + conj(dperk1)*conj(perk2)*conj(dperk3)
    per2   = (-k1^2-k2^2-k3^2)*conj(perk1)*conj(perk2)*conj(perk3)
    per3   = (((k2^2-k1^2-k3^2)/k1^2+(k1^2-k2^2-k3^2)/k2^2)*conj(dperk1)*conj(dperk2)*conj(perk3)
	           +((k2^2-k3^2-k1^2)/k3^2+(k3^2-k2^2-k1^2)/k2^2)*conj(dperk3)*conj(dperk2)*conj(perk1)
	           +((k3^2-k1^2-k2^2)/k1^2+(k1^2-k3^2-k2^2)/k3^2)*conj(dperk1)*conj(dperk3)*conj(perk2))
    per4   = conj(perk1)*conj(perk2)*conj(dperk3) + conj(dperk1)*conj(perk2)*conj(perk3) + conj(perk1)*conj(dperk2)*conj(perk3)
    per5   = (((k2^2-k1^2-k3^2)/k1^2+(k1^2-k2^2-k3^2)/k2^2)*conj(dperk1)*conj(dperk2)*conj(perk3)
	           +((k2^2-k3^2-k1^2)/k3^2+(k3^2-k2^2-k1^2)/k2^2)*conj(dperk3)*conj(dperk2)*conj(perk1)
	           +((k3^2-k1^2-k2^2)/k1^2+(k1^2-k3^2-k2^2)/k3^2)*conj(dperk1)*conj(dperk3)*conj(perk2))
    per6   = (k3^2*(k3^2-k1^2-k2^2)/k1^2/k2^2*conj(dperk1)*conj(dperk2)*conj(perk3)
                +k1^2*(k1^2-k3^2-k2^2)/k3^2/k2^2*conj(dperk3)*conj(dperk2)*conj(perk1)
                +k2^2*(k2^2-k1^2-k3^2)/k1^2/k3^2*conj(dperk1)*conj(dperk3)*conj(perk2))
    return sum([-4BG*per1,2BG2*per2,2BG3*per3,-2BG4*per4,-BG5*per5/2,-BG6*per6/2])
end
function UnIntτ(a ,H ,ϵH ,η,k1,k2,k3,Reζτ1 ,Reζτ2 ,Reζτ3 ,Imζτ1 ,Imζτ2 ,Imζτ3 ,Redζτ1 ,Redζτ2 ,Redζτ3 ,Imdζτ1 ,Imdζτ2 ,Imdζτ3 )
    perk1  = complex(Reζτ1,Imζτ1)
    perk2  = complex(Reζτ2,Imζτ2)
    perk3  = complex(Reζτ3,Imζτ3)
    dperk1 = complex(Redζτ1,Imdζτ1)
    dperk2 = complex(Redζτ2,Imdζτ2)
    dperk3 = complex(Redζτ3,Imdζτ3)
    term7  = (2a^2*ϵH*η*(conj(perk1)*conj(perk2)*conj(dperk3)
		      +conj(dperk1)*conj(perk2)*conj(perk3)+conj(perk1)*conj(dperk2)*conj(perk3)))
    term8  = (2*(a/H*conj(perk1)*conj(perk2)*conj(perk3))*(54*(a*H)^2+2*(1-ϵH)*(-(k1^2+k2^2+k3^2)/2)
		      +1/2/(a*H)^2*(k3^2*(k3^2-k1^2-k2^2)/2+k2^2*(k2^2-k1^2-k3^2)/2+k1^2*(k1^2-k2^2-k3^2)/2)))
    term9  = (-2*(ϵH/2/H^2*conj(perk1)*conj(perk2)*conj(dperk3)
        		*(k1^2+k2^2-(k2^2-k1^2-k3^2)^2/4/k3^2-(k1^2-k2^2-k3^2)^2/4/k3^2))
        		-2*(ϵH/2/H^2*conj(dperk1)*conj(perk2)*conj(perk3)
        		*(k3^2+k2^2-(k2^2-k3^2-k1^2)^2/4/k1^2-(k3^2-k2^2-k1^2)^2/4/k1^2))
        		-2*(ϵH/2/H^2*conj(perk1)*conj(dperk2)*conj(perk3)
        		*(k1^2+k3^2-(k3^2-k1^2-k2^2)^2/4/k2^2-(k1^2-k3^2-k2^2)^2/4/k2^2)))
    term10 = (2*(a*ϵH/H*conj(perk1)*conj(dperk2)*conj(dperk3)*(2-ϵH+ϵH*(k1^2-k2^2-k3^2)^2/4/k2^2/k3^2))
        		+2*(a*ϵH/H*conj(dperk1)*conj(perk2)*conj(dperk3)*(2-ϵH+ϵH*(k2^2-k1^2-k3^2)^2/4/k1^2/k3^2))
        		+2*(a*ϵH/H*conj(dperk1)*conj(dperk2)*conj(perk3)*(2-ϵH+ϵH*(k3^2-k2^2-k1^2)^2/4/k2^2/k1^2)))
    return [term7,term8,term9,term10]
end
function UnIntN(a,H,eps1,eps2,k1,k2,k3,rlk1,rlk2,rlk3,imk1,imk2,imk3,rlpk1,rlpk2,rlpk3,impk1,impk2,impk3)
	perk1  = complex(rlk1,imk1)
	perk2  = complex(rlk2,imk2)
	perk3  = complex(rlk3,imk3)
	dperk1 = complex(rlpk1,impk1)
	dperk2 = complex(rlpk2,impk2)
	dperk3 = complex(rlpk3,impk3)
	term7=2*a^3*H*eps1*eps2*(conj(perk1)*conj(perk2)*conj(dperk3)
		+conj(dperk1)*conj(perk2)*conj(perk3)+conj(perk1)*conj(dperk2)*conj(perk3))
	term8=2*(a/H*conj(perk1)*conj(perk2)*conj(perk3))*(54*(a*H)^2+2*(1-eps1)*(-(k1^2+k2^2+k3^2)/2)
		+1/2/(a*H)^2*(k3^2*(k3^2-k1^2-k2^2)/2+k2^2*(k2^2-k1^2-k3^2)/2+k1^2*(k1^2-k2^2-k3^2)/2))
	term9=(-2*(a*eps1/2/H*conj(perk1)*conj(perk2)*conj(dperk3)
		*(k1^2+k2^2-(k2^2-k1^2-k3^2)^2/4/k3^2-(k1^2-k2^2-k3^2)^2/4/k3^2))
		-2*(a*eps1/2/H*conj(dperk1)*conj(perk2)*conj(perk3)
		*(k3^2+k2^2-(k2^2-k3^2-k1^2)^2/4/k1^2-(k3^2-k2^2-k1^2)^2/4/k1^2))
		-2*(a*eps1/2/H*conj(perk1)*conj(dperk2)*conj(perk3)
		*(k1^2+k3^2-(k3^2-k1^2-k2^2)^2/4/k2^2-(k1^2-k3^2-k2^2)^2/4/k2^2)))
	term10=2*(a^3*eps1*H*conj(perk1)*conj(dperk2)*conj(dperk3)*(2-eps1+eps1*(k1^2-k2^2-k3^2)^2/4/k2^2/k3^2))
		+2*(a^3*eps1*H*conj(dperk1)*conj(perk2)*conj(dperk3)*(2-eps1+eps1*(k2^2-k1^2-k3^2)^2/4/k1^2/k3^2))
		+2*(a^3*eps1*H*conj(dperk1)*conj(dperk2)*conj(perk3)*(2-eps1+eps1*(k3^2-k2^2-k1^2)^2/4/k2^2/k1^2))
	return [term7,0*term8,term9,term10]
end
function FNLIntegrate(k1,k2,k3,y::Array,ae,ai)
    kmin  = min(k1,k2,k3)
    kmax  = max(k1,k2,k3)
    ksq   = sort([k1,k2,k3])
    # temp1 = EvolveK(y,k1,kmax,PcsPert,ae,ai))
    # temp2 = EvolveK(y,k2,kmax,PcsPert,ae,ai))
    # temp3 = EvolveK(y,k3,kmax,PcsPert,ae,ai))

    if k1==k2==k3
        temp1  = EvolveK(y,k1,kmax,PcsPert,ae,ai)
        τ1     = temp1[6]
        Reζτ1  = itp.UnivariateSpline(τ1,temp1[2][1],s=0)#interpolate(τ1,temp1[2][1],Gridded(Linear()))
        Imζτ1  = itp.UnivariateSpline(τ1,temp1[2][2],s=0)
        Redζτ1 = itp.UnivariateSpline(τ1,temp1[2][3],s=0)
        Imdζτ1 = itp.UnivariateSpline(τ1,temp1[2][4],s=0)
        aτ     = itp.UnivariateSpline(τ1,temp1[4][1],s=0)
        Hτ     = itp.UnivariateSpline(τ1,temp1[4][2],s=0)
        ϵHτ    = itp.UnivariateSpline(τ1,temp1[4][3],s=0)
        ητ     = itp.UnivariateSpline(τ1,temp1[4][4],s=0)
        dητ    = itp.UnivariateSpline(τ1,temp1[4][5],s=0)
        N1     = temp1[1]
        ReζN1  = itp.UnivariateSpline(N1,temp1[3][1],s=0)
        ImζN1  = itp.UnivariateSpline(N1,temp1[3][2],s=0)
        RedζN1 = itp.UnivariateSpline(N1,temp1[3][3],s=0)
        ImdζN1 = itp.UnivariateSpline(N1,temp1[3][4],s=0)
        aN     = itp.UnivariateSpline(N1,temp1[5][1],s=0)
        HN     = itp.UnivariateSpline(N1,temp1[5][2],s=0)
        ϵHN    = itp.UnivariateSpline(N1,temp1[5][3],s=0)
        ηN     = itp.UnivariateSpline(N1,temp1[5][4],s=0)
        dηN    = itp.UnivariateSpline(N1,temp1[5][5],s=0)
        τ2=τ1
        Reζτ2,Imζτ2,Redζτ2,Imdζτ2=Reζτ1,Imζτ1,Redζτ1,Imdζτ1
        N2=N1
        ReζN2,ImζN2,RedζN2,ImdζN2=ReζN1,ImζN1,RedζN1,ImdζN1
        τ3=τ1
        Reζτ3,Imζτ3,Redζτ3,Imdζτ3=Reζτ1,Imζτ1,Redζτ1,Imdζτ1
        N3=N1
        ReζN3,ImζN3,RedζN3,ImdζN3=ReζN1,ImζN1,RedζN1,ImdζN1
    elseif k1!=k2!=k3 && k1!=k3
        temp3  = EvolveK(y,ksq[3],kmax,PcsPert,ae,ai)
        τ3     = temp3[6]
        Reζτ3  = itp.UnivariateSpline(τ3,temp3[2][1],s=0)
        Imζτ3  = itp.UnivariateSpline(τ3,temp3[2][2],s=0)
        Redζτ3 = itp.UnivariateSpline(τ3,temp3[2][3],s=0)
        Imdζτ3 = itp.UnivariateSpline(τ3,temp3[2][4],s=0)
        aτ     = itp.UnivariateSpline(τ3,temp3[4][1],s=0)
        Hτ     = itp.UnivariateSpline(τ3,temp3[4][2],s=0)
        ϵHτ    = itp.UnivariateSpline(τ3,temp3[4][3],s=0)
        ητ     = itp.UnivariateSpline(τ3,temp3[4][4],s=0)
        dητ    = itp.UnivariateSpline(τ3,temp3[4][5],s=0)
        N3     = temp3[1]
        ReζN3  = itp.UnivariateSpline(N3,temp3[3][1],s=0)
        ImζN3  = itp.UnivariateSpline(N3,temp3[3][2],s=0)
        RedζN3 = itp.UnivariateSpline(N3,temp3[3][3],s=0)
        ImdζN3 = itp.UnivariateSpline(N3,temp3[3][4],s=0)
        aN     = itp.UnivariateSpline(N3,temp3[5][1],s=0)
        HN     = itp.UnivariateSpline(N3,temp3[5][2],s=0)
        ϵHN    = itp.UnivariateSpline(N3,temp3[5][3],s=0)
        ηN     = itp.UnivariateSpline(N3,temp3[5][4],s=0)
        dηN    = itp.UnivariateSpline(N3,temp3[5][5],s=0)

        temp2  = EvolveK(y,ksq[2],kmax,PcsPert,ae,ai)
        τ2     = temp2[6]
        Reζτ2  = itp.UnivariateSpline(τ2,temp2[2][1],s=0)
        Imζτ2  = itp.UnivariateSpline(τ2,temp2[2][2],s=0)
        Redζτ2 = itp.UnivariateSpline(τ2,temp2[2][3],s=0)
        Imdζτ2 = itp.UnivariateSpline(τ2,temp2[2][4],s=0)
        N2     = temp2[1]
        ReζN2  = itp.UnivariateSpline(N2,temp2[3][1],s=0)
        ImζN2  = itp.UnivariateSpline(N2,temp2[3][2],s=0)
        RedζN2 = itp.UnivariateSpline(N2,temp2[3][3],s=0)
        ImdζN2 = itp.UnivariateSpline(N2,temp2[3][4],s=0)

        temp1  = EvolveK(y,ksq[1],kmax,PcsPert,ae,ai)
        τ1     = temp1[6]
        Reζτ1  = itp.UnivariateSpline(τ1,temp1[2][1],s=0)
        Imζτ1  = itp.UnivariateSpline(τ1,temp1[2][2],s=0)
        Redζτ1 = itp.UnivariateSpline(τ1,temp1[2][3],s=0)
        Imdζτ1 = itp.UnivariateSpline(τ1,temp1[2][4],s=0)
        N1     = temp1[1]
        ReζN1  = itp.UnivariateSpline(N1,temp1[3][1],s=0)
        ImζN1  = itp.UnivariateSpline(N1,temp1[3][2],s=0)
        RedζN1 = itp.UnivariateSpline(N1,temp1[3][3],s=0)
        ImdζN1 = itp.UnivariateSpline(N1,temp1[3][4],s=0)
    else
        temp1  = EvolveK(y,ksq[1],kmax,PcsPert,ae,ai)
        τ1     = temp1[6]
        Reζτ1  = itp.UnivariateSpline(τ1,temp1[2][1],s=0)
        Imζτ1  = itp.UnivariateSpline(τ1,temp1[2][2],s=0)
        Redζτ1 = itp.UnivariateSpline(τ1,temp1[2][3],s=0)
        Imdζτ1 = itp.UnivariateSpline(τ1,temp1[2][4],s=0)
        N1     = temp1[1]
        ReζN1  = itp.UnivariateSpline(N1,temp1[3][1],s=0)
        ImζN1  = itp.UnivariateSpline(N1,temp1[3][2],s=0)
        RedζN1 = itp.UnivariateSpline(N1,temp1[3][3],s=0)
        ImdζN1 = itp.UnivariateSpline(N1,temp1[3][4],s=0)

        temp3  = EvolveK(y,ksq[3],kmax,PcsPert,ae,ai)
        τ3     = temp3[6]
        Reζτ3  = itp.UnivariateSpline(τ3,temp3[2][1],s=0)
        Imζτ3  = itp.UnivariateSpline(τ3,temp3[2][2],s=0)
        Redζτ3 = itp.UnivariateSpline(τ3,temp3[2][3],s=0)
        Imdζτ3 = itp.UnivariateSpline(τ3,temp3[2][4],s=0)
        aτ     = itp.UnivariateSpline(τ3,temp3[4][1],s=0)
        Hτ     = itp.UnivariateSpline(τ3,temp3[4][2],s=0)
        ϵHτ    = itp.UnivariateSpline(τ3,temp3[4][3],s=0)
        ητ     = itp.UnivariateSpline(τ3,temp3[4][4],s=0)
        dητ    = itp.UnivariateSpline(τ3,temp3[4][5],s=0)
        N3     = temp3[1]
        ReζN3  = itp.UnivariateSpline(N3,temp3[3][1],s=0)
        ImζN3  = itp.UnivariateSpline(N3,temp3[3][2],s=0)
        RedζN3 = itp.UnivariateSpline(N3,temp3[3][3],s=0)
        ImdζN3 = itp.UnivariateSpline(N3,temp3[3][4],s=0)
        aN     = itp.UnivariateSpline(N3,temp3[5][1],s=0)
        HN     = itp.UnivariateSpline(N3,temp3[5][2],s=0)
        ϵHN    = itp.UnivariateSpline(N3,temp3[5][3],s=0)
        ηN     = itp.UnivariateSpline(N3,temp3[5][4],s=0)
        dηN    = itp.UnivariateSpline(N3,temp3[5][5],s=0)
        if ksq[2]==ksq[1]
            τ2=τ1
            Reζτ2,Imζτ2,Redζτ2,Imdζτ2=Reζτ1,Imζτ1,Redζτ1,Imdζτ1
            N2=N1
            ReζN2,ImζN2,RedζN2,ImdζN2=ReζN1,ImζN1,RedζN1,ImdζN1
        else
            τ2=τ3
            Reζτ2,Imζτ2,Redζτ2,Imdζτ2=Reζτ3,Imζτ3,Redζτ3,Imdζτ3
            N2=N3
            ReζN2,ImζN2,RedζN2,ImdζN2=ReζN3,ImζN3,RedζN3,ImdζN3
        end
    end
    τstop  = last(τ3)
    τstart = τ3[1]
    Nstop  = last(N3)
    Nstart = N3[1]
    Rζ1,Rζ2,Rζ3=ReζN1(Nstop)[1],ReζN2(Nstop)[1],ReζN3(Nstop)[1]
    Iζ1,Iζ2,Iζ3=ImζN1(Nstop)[1],ImζN2(Nstop)[1],ImζN3(Nstop)[1]
    λreal=(Rζ2*Rζ3*Iζ1+Rζ1*Rζ3*Iζ2+Rζ1*Rζ2*Iζ3-Iζ1*Iζ2*Iζ3)
    λimag=(Rζ1*Rζ2*Rζ3-Rζ3*Iζ1*Iζ2-Rζ2*Iζ1*Iζ3-Rζ1*Iζ2*Iζ3)
    ζ1,ζ2,ζ3=complex(Rζ1,Iζ1),complex(Rζ2,Iζ2),complex(Rζ3,Iζ3)
    𝒫1,𝒫2,𝒫3=abs(ζ1)^2,abs(ζ2)^2,abs(ζ3)^2
    𝒫12,𝒫13,𝒫23=𝒫1*𝒫2,𝒫1*𝒫3,𝒫2*𝒫3

    # display("τstart=$τstart,Nstop=$Nstop")
    # display("$([aN(Nstop)[1],HN(Nstop)[1],ϵHN(Nstop)[1],ηN(Nstop)[1],ReζN1(Nstop)[1],ReζN2(Nstop)[1],ReζN3(Nstop)[1],ImζN1(Nstop)[1]
    #                 ,ImζN2(Nstop)[1],ImζN3(Nstop)[1],RedζN1(Nstop)[1],RedζN2(Nstop)[1],RedζN3(Nstop)[1],ImdζN1(Nstop)[1],ImdζN2(Nstop)[1],ImdζN3(Nstop)[1]])")
    # display(aτ(τstart)[1])
    function IntKernelτ(x)
        temp=Kernelτ(aτ(x)[1],Hτ(x)[1],ϵHτ(x)[1],dητ(x)[1],k1,k2,k3,Reζτ1(x)[1],Reζτ2(x)[1],Reζτ3(x)[1],
        Imζτ1(x)[1],Imζτ2(x)[1],Imζτ3(x)[1],Redζτ1(x)[1],Redζτ2(x)[1],Redζτ3(x)[1],Imdζτ1(x)[1],Imdζτ2(x)[1],Imdζτ3(x)[1])
        return exp(κ*kmax*(x-τstop))*(λreal*real(temp)+λimag*imag(temp))
    end
    function IntKernelN(x)
        temp=KernelN(aN(x)[1],HN(x)[1],ϵHN(x)[1],dηN(x)[1],k1,k2,k3,ReζN1(x)[1],ReζN2(x)[1],ReζN3(x)[1],
        ImζN1(x)[1],ImζN2(x)[1],ImζN3(x)[1],RedζN1(x)[1],RedζN2(x)[1],RedζN3(x)[1],ImdζN1(x)[1],ImdζN2(x)[1],ImdζN3(x)[1])
        return λreal*real(temp)+λimag*imag(temp)
    end
    # UnIntτTemp = UnIntτ(aτ(τstart)[1],Hτ(τstart)[1],ϵHτ(τstart)[1],ητ(τstart)[1],k1,k2,k3,Reζτ1(τstart)[1],Reζτ2(τstart)[1],Reζτ3(τstart)[1],Imζτ1(τstart)[1]
    #                 ,Imζτ2(τstart)[1],Imζτ3(τstart)[1],Redζτ1(τstart)[1],Redζτ2(τstart)[1],Redζτ3(τstart)[1],Imdζτ1(τstart)[1],Imdζτ2(τstart)[1],Imdζτ3(τstart)[1])
    UnIntNTemp = UnIntN(aN(Nstop)[1],HN(Nstop)[1],ϵHN(Nstop)[1],ηN(Nstop)[1],k1,k2,k3,ReζN1(Nstop)[1],ReζN2(Nstop)[1],ReζN3(Nstop)[1],ImζN1(Nstop)[1]
                    ,ImζN2(Nstop)[1],ImζN3(Nstop)[1],RedζN1(Nstop)[1],RedζN2(Nstop)[1],RedζN3(Nstop)[1],ImdζN1(Nstop)[1],ImdζN2(Nstop)[1],ImdζN3(Nstop)[1])
    UnInt      = sum(UnIntNTemp)
    UnIntTerm = λreal*real(UnInt)+λimag*imag(UnInt)
    abstol=1e-4*max(𝒫12,𝒫13,𝒫23,UnIntTerm)
    ifprint==1 ? display("UnInt=$UnInt,abstol=$abstol,UnIntTerm=$UnIntTerm") : nothing
    # ReIntτ,=hquadrature(ReIntKernelτ,τstart,τstop,reltol=1e-15,abstol=1e-22)
    # ImIntτ,=hquadrature(ImIntKernelτ,τstart,τstop,reltol=1e-15,abstol=1e-22)
    # ReIntN,=hquadrature(ReIntKernelN,Nstart,Nstop,reltol=1e-15,abstol=1e-22)
    # ImIntN,=hquadrature(ImIntKernelN,Nstart,Nstop,reltol=1e-15,abstol=1e-22)
    # Intτ=complex(ReIntτ,ImIntτ)
    # IntN=complex(ReIntN,ImIntN)
    # display("UnInt=$UnInt")
    Intτ, = quadgk(IntKernelτ,τstart,τstop,rtol=1e-16,atol=abstol)
    IntN, = quadgk(IntKernelN,Nstart,Nstop,rtol=1e-16,atol=abstol)
    # if 1e12<=kmax<=10e12
    #     Intτ, = quadgk(IntKernelτ,τstart,τstop,rtol=1e-16,atol=1e-77)
    #     IntN, = quadgk(IntKernelN,Nstart,Nstop,rtol=1e-16,atol=1e-64)
    # elseif 10e12<kmax<=1000e12
    #     Intτ, = quadgk(IntKernelτ,τstart,τstop,rtol=1e-16,atol=1e-26)
    #     IntN, = quadgk(IntKernelN,Nstart,Nstop,rtol=1e-16,atol=1e-20)
    # else
    #     Intτ, = quadgk(IntKernelτ,τstart,τstop)
    #     IntN, = quadgk(IntKernelN,Nstart,Nstop)
    # end
        # Intτ=itg.romberg(IntKernelτ,τstart,τstop,rtol=1e-16,tol=1e-20,divmax=16)
        # IntN=itg.romberg(IntKernelN,Nstart,Nstop,rtol=1e-16,tol=1e-16,divmax=16)
        # Intτ, = quadgk(IntKernelτ,τstart,τstop)
    # ifprint==1 ? display("UnInt=$UnInt,Intτ=0,IntN=$IntN") : nothing
    ifprint==1 ? display("UnIntTerm=$UnIntTerm,Intτ=$Intτ,IntN=$IntN") : nothing
    Bi=Intτ+IntN+UnIntTerm#(UnInt+IntN))
    fnl=2Bi/(𝒫12+𝒫13+𝒫23)
    #display("fnl=$fnl,Ps1=$𝒫1,Ps2=$𝒫2,Ps3=$𝒫3")
    return [fnl,Bi,𝒫1,ksq[1]]
end
function FNLProcedure(k1,k2,k3)
    kmin  = min(k1,k2,k3)
    kmax  = max(k1,k2,k3)
    maxtime = 1000
    Dϕpivot,Hpivot=Attractor(ϕpivot,PcsBack,PcsPivot)#Dϕ=dϕ/dt while dϕ=dϕ/dτ
    apivot=kpivot/Hpivot
    y1=[apivot,ϕpivot,apivot*Dϕpivot]
    # display(y1)
    if kmin/RatioMin>=kpivot
        y = y1
    else
        aHini=kmin/RatioMin
        atry =apivot
        ϕtry =ϕpivot
        Htry =Hpivot
        Dϕtry=Dϕpivot
        counter=0
        while atry*Htry>=aHini
            counter +=1
            if counter>=maxtime
                throw(ErrorException("When searching for an initial value of ϕ, the code couldn't converge after $maxtime iterations.
                The potential might not allow enough e-folds before reaching pivot scale."))
            end
            V,dV=Potential(ϕtry)
            g,dg=Gfun(ϕtry)
            Dϕ=Dϕtry
            H=Hubble(Dϕ,g,V)#sqrt((dϕ^2/2+V-g*dϕ^2)/3)
            ϕtry += -1.01*log(atry*Htry/aHini)*Dϕ/H
            Dϕtry,Htry=Attractor(ϕtry,PcsBack,PcsIni)
            y = [1,ϕtry,Dϕtry]
            y = EvolveBackground(y,ϕpivot,PcsBack)
            atry=apivot/y[1]
        end
        y = [atry,ϕtry,atry*Dϕtry]
    end
    dy   = Derivs(y,0,3)
    Htry = dy[1]/y[1]^2
    ifprint==1 ? display("Initially, we have a=$(y[1]),H=$Htry,ϕ=$(y[2]),k=$(Htry*y[1])") : nothing
    ai=y1[1]
    temp=ReachInf(y,PcsBack,ai)
    #display(temp)
    ae=temp[3]
    return FNLIntegrate(k1,k2,k3,y,ae,ai)
end
function FNL(k1,k2,k3)
    if k1+k2>k3 && k1+k3>k2 && k2+k3>k1
        return FNLProcedure(k1,k2,k3)
    else
        throw(ErrorException("kᵢ couldn't form a triangle."))
    end
end
function SampleK()
    k=10 .^(range(log10(5e13),stop=log10(5e14),length=200))
    return k
end
function ComputingFNL()
    fnltrack=fill([0.0],length(ktrack))
    for i in 1:length(ktrack)
        display(i)
        temp=ktrack[i]
        fnltrack[i] = FNL(temp,temp,temp)
    end
    return fnltrack
end
function Spectrum()
    maxtime = 1000
    Dϕpivot,Hpivot=Attractor(ϕpivot,PcsBack,PcsPivot)#Dϕ=dϕ/dt while dϕ=dϕ/dτ
    apivot=kpivot/Hpivot
    y1=[apivot,ϕpivot,apivot*Dϕpivot]
    aHini=kminimum/RatioMin
    atry =apivot
    ϕtry =ϕpivot
    Htry =Hpivot
    Dϕtry=Dϕpivot
    counter=0
    while atry*Htry>=aHini
        counter +=1
        if counter>=maxtime
            throw(ErrorException("When searching for an initial value of ϕ, the code couldn't converge after $maxtime iterations.
            The potential might not allow enough e-folds before reaching pivot scale."))
        end
        V,dV=Potential(ϕtry)
        g,dg=Gfun(ϕtry)
        Dϕ=Dϕtry
        H=Hubble(Dϕ,g,V)#sqrt((dϕ^2/2+V-g*dϕ^2)/3)
        ϕtry += -1.01*log(atry*Htry/aHini)*Dϕ/H
        Dϕtry,Htry=Attractor(ϕtry,PcsBack,PcsIni)
        y = [1,ϕtry,Dϕtry]
        y = EvolveBackground(y,ϕpivot,PcsBack)
        atry=apivot/y[1]
    end
    y    = [atry,ϕtry,atry*Dϕtry]
    dy   = Derivs(y,0,3)
    Htry = dy[1]/y[1]^2
    ai   = y1[1]
    ifprint==1 ? display("Initially, we have a=$(y[1]),H=$Htry,ϕ=$(y[2]),k=$(Htry*y[1])") : nothing
    ϕe,ke,ae,Ne=ReachInf(y,PcsBack,ai)
    return ComputingSpectrum(y,ae)
end
function Mode(yi::Array,k,accuracy,ae)
    y    = ReachAH(yi,k/RatioMin,PcsBack)
    dy   = Derivs(y,k,3)
    dV   = Potential(y[2])[2]
    H    = dy[1]/y[1]^2
    aH   = H*y[1]
    dϕ   = dy[2]/y[1]
    para = CheckSR(y[2],dϕ)
    ϵH,η = para
    g,   = Gfun(y[2])
    z    = y[1]*dϕ/H*sqrt(1-2g)
    zpz  = y[1]*H*(1+η/2)
    yini = [1/sqrt(2k)/z,0,-zpz/sqrt(2k)/z,-k/sqrt(2k)/z]
    y    = vcat(y,yini)
    dτ   = accuracy*2π/k
    ratio= k/aH
    while ratio>=RatioMax && y[1]<=ae
        y  = Rk4(y,11,k,dτ)
        dy = Derivs(y,k,11)
        aH = dy[1]/y[1]
        ratio = k/aH
        dτ = accuracy*2π/max(sqrt(abs(dy[6]/y[4])),k,aH,abs(dy[3]/y[3]))
    end
    ζ=y[4]^2+y[5]^2
    return k^3/(2π^2)*ζ
end
function ComputingSpectrum(yi::Array,aend)
    y    = yi
    V,   = Potential(y[2])
    g,   = Gfun(y[2])
    dϕ   = y[3]/y[1]
    aH   = y[1]*Hubble(dϕ,g,V)
    temp = 0
    curvature = Float64[]
    if aH>=kminimum/RatioMin
        throw(ErrorException("At initial time, kmin is already outside the horizon"))
    end
    ktrack=10 .^(range(12,stop=14,length=500))
    #vcat(10 .^(range(-4,stop=12,length=200)),10 .^(range(12,stop=14,length=500)),10 .^(range(14,stop=18,length=200)))
    for i in 1:length(ktrack)
        display(i)
        y    = ReachAH(y,ktrack[i]/100,PcsBack)
        temp = Mode(y,ktrack[i],PcsBack,aend)
        push!(curvature,temp)
    end
    spectra=[ktrack,curvature]
    writedlm("Spectrum.csv",spectra, ",")
    plot(spectra)
    return spectra
end
const ktrack=SampleK()
# # const ktrack= 10 .^(range(0,stop=15,length=200))
@time temp=ComputingFNL()
temp2=[temp[i][1] for i in 1:length(ktrack)]
# # writedlm("higgs_fnl_eq.dat", [ktrack,temp2], ";")
# # plot(ktrack,temp2,xscale=:log10)
# # writedlm("fnl1.csv", temp, ";")
# # for i in 1:length(ktrack)
# #     println(ktrack[i])
# # end
# # temp=ktrack[1]
# const ktemp=1e14
# @time display(FNL(ktemp,ktemp,ktemp))
# println(Potential(ϕpivot))
# km1,km2=1e8,1e11
# nk1,nk2,nk3=10,30,60
# rk1=(km1/1)^(1/nk1)
# rk2=(km2/km1)^(1/nk2)
# rk3=(1e15/km2)^(1/nk3)
# ktrack=[1e11]
# kk=1e11
# for i in 1:(nk1+nk2+nk3)
#     if 1<=i<nk1
#         kk *=rk1
#         push!(ktrack,kk)
#     elseif nk1<=i<nk1+nk2
#         kk *=rk2
#         push!(ktrack,kk)
#     else
#         kk *= rk3
#         push!(ktrack,kk)
#     end
# end
# @time temp=Spectrum()
# ktrack=Float64[]
# temp2 = Float64[]
# for i in 1:(length(temp[1])-1)
#     a = (log(temp[2][i+1])-log(temp[2][i]))/(log(temp[1][i+1])-log(temp[1][i]))
#     push!(temp2,a)
#     push!(ktrack,temp[1][i])
# end
