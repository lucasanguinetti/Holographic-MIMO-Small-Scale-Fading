function [var_lm] = function_computeInt(l,m,Lx,Ly,k_phi1,k_phi2,k_phi3,k_phi4)

%Function Inputs
%l,m:                           discrete wavenumber frequency over (kx,ky)
%Lx,Ly:                         normalized array size
%k_phi1,k_phi2,k_phi3,k_phi4:   angular integration regions 

%Function Outputs
%var_lm:      variance of (l,m)-th Fourier coefficient in the first
%quadrant (the other are obtained by symmetry)

if m==0 && l==0
    %1st integral
    int1 = 0; %k_phi1=k_phi_2=0
    %2nd integral
    %2.1
    int21 = k_phi3;
    %2.2
    a2 = Lx;
    if a2==1
        int22 = 0;
    else
        int22 = cosInt(a2,min(k_phi3,acos(1/a2)))-cosInt(a2,0);
    end
    %sum
    int2 = int21-int22;
    %3rd integral
    %3.1
    int31 = pi/2-k_phi3;
    %3.2
    a3 = Ly;
    if a3==1
        int32 = 0;
    else
        int32 = sinInt(a3,pi/2)-sinInt(a3,max(k_phi3,asin(1/a3))); %k_phi4=pi/2
    end
    %sum
    int3 = int31-int32;
    
elseif m==0 && l~=0
    %1st integral
    int1 = 0; %k_phi1=k_phi2=0
    %2nd integral
    a21 = Lx/l;
    int21 = cosInt(a21,min(k_phi3,acos(1/a21)))-cosInt(a21,min(k_phi2,acos(1/a21)));
    %2.2
    a22 = Lx/(l+1);
    if a22==1
        int22 = 0;
    else
        int22 = cosInt(a22,min(k_phi3,acos(1/a22)))-cosInt(a22,min(k_phi2,acos(1/a22)));
    end
    %sum
    int2 = int21-int22;
    %3rd integral
    %3.1
    a31 = Lx/l;
    int31 = cosInt(a31,min(k_phi4,acos(1/a31)))-cosInt(a31,min(k_phi3,acos(1/a31)));
    %3.2
    a32 = Ly/(m+1);
    if a32==1
        int32 = 0;
    else
        int32 = sinInt(a32,max(k_phi4,asin(1/a32)))-sinInt(a32,max(k_phi3,asin(1/a32)));
    end
    %sum
    int3 = int31-int32;
    
elseif m~=0 && l==0
    %1st integral
    %1.1
    a11 = Ly/m;
    int11 = sinInt(a11,max(k_phi2,asin(1/a11)))-sinInt(a11,max(k_phi1,asin(1/a11)));
    %1.2
    a12 = Lx/(l+1);
    if a12==1
        int12 = 0;
    else
        int12 = cosInt(a12,min(k_phi2,acos(1/a12)))-cosInt(a12,min(k_phi1,acos(1/a12)));
    end
    %sum
    int1 = int11-int12;
    %2nd integral
    %2.1
    a21 = Ly/m;
    int21 = sinInt(a21,max(k_phi3,asin(1/a21)))-sinInt(a21,max(k_phi2,asin(1/a21)));
    %2.2
    a22 = Ly/(m+1);
    if a22==1
        int22 = 0;
    else
        int22 = sinInt(a22,max(k_phi3,asin(1/a22)))-sinInt(a22,max(k_phi2,asin(1/a22)));
    end
    %sum
    int2 = int21-int22;
    %3rd integral
    int3 = 0; %k_phi3=k_phi4=pi/2
    
else
    %1st integral
    %1.1
    a11 = Ly/m;
    int11 = sinInt(a11,max(k_phi2,asin(1/a11)))-sinInt(a11,max(k_phi1,asin(1/a11)));
    %1.2
    a12 = Lx/(l+1);
    if a12==1
        int12 = 0;
    else
        int12 = cosInt(a12,min(k_phi2,acos(1/a12)))-cosInt(a12,min(k_phi1,acos(1/a12)));
    end
    %sum
    int1 = int11-int12;
    %2nd integral
    if l>=m
        %2.1
        a21 = Lx/l;
        int21 = cosInt(a21,min(k_phi3,acos(1/a21)))-cosInt(a21,min(k_phi2,acos(1/a21)));
        %2.2
        a22 = Lx/(l+1);
        if a22==1
            int22 = 0;
        else
            int22 = cosInt(a22,min(k_phi3,acos(1/a22)))-cosInt(a22,min(k_phi2,acos(1/a22)));
        end
    else
        %2.1
        a21 = Ly/m;
        int21 = sinInt(a21,max(k_phi3,asin(1/a21)))-sinInt(a21,max(k_phi2,asin(1/a21)));
        %2.2
        a22 = Ly/(m+1);
        if a22==1
            int22 = 0;
        else
            int22 = sinInt(a22,max(k_phi3,asin(1/a22)))-sinInt(a22,max(k_phi2,asin(1/a22)));
        end
    end
    %sum
    int2 = int21-int22;
    %3rd integral
    %3.1
    a31 = Lx/l;
    int31 = cosInt(a31,min(k_phi4,acos(1/a31)))-cosInt(a31,min(k_phi3,acos(1/a31)));
    %3.2
    a32 = Ly/(m+1);
    if a32==1
        int32 = 0;
    else
        int32 = sinInt(a32,max(k_phi4,asin(1/a32)))-sinInt(a32,max(k_phi3,asin(1/a32)));
    end
    %sum
    int3 = int31-int32;
end

%sum up all contributions
var_lm = real(int1+int2+int3);

end
