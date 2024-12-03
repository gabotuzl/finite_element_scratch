close
clc
clear

% TO USE THIS SCRIPT, FOLLOW THE SAME INDICATIONS PROVIDED FOR THE PYTHON SCRIPTS IN THE README FILE. THE ONLY DIFFERENCE IS THAT ALL THE FUNCTIONS ARE DEFINED IN THIS ONE SCRIPT FOR THE MATLAB CASE.


syms x chi

%Parameters
E=600e6; 
A=0.00589;
I=0.00004545;

%Node coordinates
NodePositions=[0 0; 1 0; 2 0; 3 0; 4 0; 5 0; 6 0; 7 0; 8 0; 0.5 0.5; 3.5 0.5; 4 0.5; 4.5 0.5; 7.5 0.5; 1 1; 3 1; 4 1; 5 1; 7 1; 1.5 1.5; 2.5 1.5;...
    4 1.5; 5.5 1.5; 6.5 1.5; 2 2; 3 2; 4 2; 5 2; 6 2];

%Connectivity matrix
connectivity=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 1 10; 11 5; 12 5; 5 13; 14 9; 10 15; 16 11; 17 12; 13 18; 19 14; 15 20; 21 16;...
    22 17; 18 23; 24 19; 20 25; 25 21; 27 22; 23 29; 29 24; 25 26; 26 27; 27 28; 28 29];

%Create element positions
for i=1:length(connectivity)
    a=connectivity(i,1);
    b=connectivity(i,2);
    elements(i,:)=[NodePositions(a,:) NodePositions(b,:)];
end

%Global matrix parameters
NodeDOF=3;
TotalDOF=NodeDOF*length(NodePositions);

%Force vector
GlobalForces=zeros(TotalDOF,1);
LoadedNode1=3;
LoadedNode2=7;
GlobalForces((LoadedNode1-1)*NodeDOF+2)=-50000;
GlobalForces((LoadedNode2-1)*NodeDOF+2)=-50000;

%Fixed nodes (1 means that that node is fixed)
FixedNodes=[1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

for k=1:length(elements)
    v1=v1Func(k,elements);
    v2=v2Func(k,elements);
    ElementAngle=(atan2((elements(k,4)-elements(k,2)),(elements(k,3)-elements(k,1))));
    he=sqrt((elements(k,1)-elements(k,3))^2+(elements(k,2)-elements(k,4))^2); %Element length
    
    %Matrix K1
    for i=1:2
        for j=1:2
            d1=diff(v1(i),chi);
            d2=diff(v1(j),chi);
            K1(i,j)=E*A*int((d1*d2),chi,0,1);
        end
    end
    
    %Matrix K2
    for i=1:4
        for j=1:4
            dd1=diff(diff(v2(i),chi),chi);
            dd2=diff(diff(v2(j),chi),chi);
            K2(i,j)=E*I*int((dd1*dd2),chi,0,1);
        end
    end

    %Construction of the local stiffness matrix KTotal
    Kinter=[K1 zeros(2,4); zeros(4,2) K2];
    KTotal=[Kinter(1,:);Kinter(3,:);Kinter(4,:);Kinter(2,:);Kinter(5,:);Kinter(6,:)];
    Kinter=KTotal;
    KTotal=[Kinter(:,1) Kinter(:,3) Kinter(:,4) Kinter(:,2) Kinter(:,5) Kinter(:,6)];
    
    %Rotating the local stiffness matrix in the global reference frame
    Rot=[cosd(ElementAngle) sin(ElementAngle) 0 0 0 0; -sin(ElementAngle) cos(ElementAngle) 0 0 0 0; 0 0 1 0 0 0; 0 0 0 cos(ElementAngle) sin(ElementAngle) 0; 0 0 0 -sin(ElementAngle) cos(ElementAngle) 0; 0 0 0 0 0 1];
    KTotal=Rot'*KTotal*Rot;
    
    %Storing the local stiffness matrix of the element
    LocalMatrices{k}=KTotal;
end   

%Construction of the global stiffness matrix GlobalMatrix
GlobalMatrix=zeros(TotalDOF,TotalDOF);
for e=1:length(connectivity)
    i=connectivity(e,1);
    j=connectivity(e,2);
    
    %Positioning at the top-left corner of the submatrix
    iG=(i-1)*NodeDOF+1;
    jG=(j-1)*NodeDOF+1;
    
    %Insert values from submatrices into the global_matrix matrix
    GlobalMatrix(iG:iG+2,iG:iG+2)=GlobalMatrix(iG:iG+2,iG:iG+2)+LocalMatrices{e}(1:3,1:3);
    GlobalMatrix(iG:iG+2,jG:jG+2)=GlobalMatrix(iG:iG+2,jG:jG+2)+LocalMatrices{e}(1:3,4:6);
    GlobalMatrix(jG:jG+2,iG:iG+2)=GlobalMatrix(jG:jG+2,iG:iG+2)+LocalMatrices{e}(4:6,1:3);
    GlobalMatrix(jG:jG+2,jG:jG+2)=GlobalMatrix(jG:jG+2,jG:jG+2)+LocalMatrices{e}(4:6,4:6);
end

%Creating new matrix variable, to accomodate for the boundary conditions
BC_GlobalMatrix=GlobalMatrix;
BC_GlobalForces=GlobalForces;

%Iterate over FixedNodes in reverse order (Evaluating the boundary conditions)
for i=1:length(FixedNodes)
    %Start from the end of FixedNodes
    k=length(FixedNodes)-i+1;
    %Deletes rows and columns associated to the fixed nodes
    if FixedNodes(k)==1
        M=(k-1)*NodeDOF+1;
        BC_GlobalMatrix(M:M+2,:)=[];
        BC_GlobalMatrix(:,M:M+2)=[];
        BC_GlobalForces(M:M+2)=[];
    end
end

%Calculate the free degrees of freedom
FreeDisplacements=BC_GlobalMatrix\BC_GlobalForces;

%Rebuild the vector of degrees of freedom, including fixed nodes
TotalDisplacements=zeros(TotalDOF,1);
j=1;
for i=1:length(FixedNodes)
    k=(i-1)*NodeDOF+1;
    if FixedNodes(i)==1
        TotalDisplacements(k)=0;
        TotalDisplacements(k+1)=0;
        TotalDisplacements(k+2)=0;
    else
        TotalDisplacements(k)=FreeDisplacements(j);
        TotalDisplacements(k+1)=FreeDisplacements(j+1);
        TotalDisplacements(k+2)=FreeDisplacements(j+2);
        j=j+3;
    end
end

%Calculating the force vector, including the fixed nodes 
NodeForces=GlobalMatrix*TotalDisplacements;
Grafica(elements,TotalDisplacements,NodeDOF,connectivity)
[MaxDeformacion,MaxEsfuerzo]=PostAnalisis(elements,TotalDisplacements,NodeDOF,E,connectivity,he)

function v1=v1Func(Num,elements)
    syms chi
    
    he=sqrt((elements(Num,1)-elements(Num,3))^2+(elements(Num,2)-elements(Num,4))^2);
    
    v1(1)=1-chi;
    v1(2)=chi;
end

function v2=v2Func(Num,elements)
    syms chi
    
    he=sqrt((elements(Num,1)-elements(Num,3))^2+(elements(Num,2)-elements(Num,4))^2);

    v2(1)=1-3*chi^2+2*chi^3;
    v2(2)=-he*chi*(1-chi)^2;
    v2(3)=3*chi^2-2*chi^3;
    v2(4)=-he*chi*(chi^2-chi);
end

function Grafica(elements,TotalDisplacements,NodeDOF,connectivity)
    figure(1)
    hold on
    size=30;
    for i=1:length(elements)
        
    %Initial node positions
    plot([elements(i,1) elements(i,3)],[elements(i,2) elements(i,4)], 'g-', 'LineWidth',1);
    scatter(elements(i,1),elements(i,2),size)
    scatter(elements(i,3),elements(i,4),size)
    %Final node positions
    c=connectivity(i,:);
    ka=(c(1)-1)*NodeDOF+1;
    kb=(c(2)-1)*NodeDOF+1;
    PosNueva=[TotalDisplacements(ka)+elements(i,1) TotalDisplacements(ka+1)+elements(i,2) TotalDisplacements(kb)+elements(i,3) TotalDisplacements(kb+1)+elements(i,4)];
    plot([PosNueva(1) PosNueva(3)],[PosNueva(2) PosNueva(4)], 'r-', 'LineWidth',1);
    scatter(PosNueva(1),PosNueva(2),size)
    scatter(PosNueva(3),PosNueva(4),size)
    end
    
end

function [MaxDeformacion,MaxEsfuerzo]=PostAnalisis(elements,TotalDisplacements,NodeDOF,E,connectivity,he)
    syms chi
    for k=1:length(elements)
        v1=v1Func(k,elements);
        v2=v2Func(k,elements);
        
        a=connectivity(k,1);
        b=connectivity(k,2);
        
        eps0=TotalDisplacements((a-1)*NodeDOF+1)*diff(v1(1),chi)+TotalDisplacements((b-1)*NodeDOF+1)*diff(v1(2),chi);
        eps1=TotalDisplacements((a-1)*NodeDOF+2)*diff(diff(v2(1),chi),chi)+TotalDisplacements((a-1)*NodeDOF+3)*diff(diff(v2(2),chi),chi)+...
            TotalDisplacements((b-1)*NodeDOF+2)*diff(diff(v2(3),chi),chi)+TotalDisplacements((b-1)*NodeDOF+3)*diff(diff(v2(4),chi),chi);
        deformacion=eps0+eps1;
        
        esfuerzo=E*eps0-E*(0.203/2)*eps1;
        
        for i=1:11
            gradDeformacion(k,i)=subs(deformacion,chi,i/10-1/10);
            gradEsfuerzo(k,i)=subs(esfuerzo,chi,i/10);
        end
        MaxDeformacion(k)=double(max(gradDeformacion(k,:)));
        MaxEsfuerzo(k)=double(max(gradEsfuerzo(k,:)))
    end
end        
        
