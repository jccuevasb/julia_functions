function [countovf, Uaver, sortu, sorty, ypperturbated]=umean_profiles_centroid(yplus,uplus,Np,gp,wvf)
%----yplus=wall normal positions without perturbate
%----Np= number of assembles Np=0, master profile
%----gp= number of grid points
%----wvf= width of vortical fissure
%----Output paramters
%----ustep=row vector of u velocities of 5481 elements
%----ystep=row vector of y+ wall normal positions of 5481 elements

%---Initialize vectors
dy=wvf;
countovf=0;                                      %--counter overlapping 
nypl=length(yplus);                              %--number of elements yplus
%ul1vf=0.5;                                      %--Velocity @ the left edge first vortical fissure default ul1vf=0.5
ul1vf=1.0;                                       %--Velocity @ the left edge first vortical fissure
ypro=zeros(nypl,1);
ypperturbated=zeros(nypl,1);
sorty=zeros(nypl,1);
sortu=zeros(nypl,1);
%gp=5481;                                        %--Grid Points
U=ones(gp-1,nypl)*NaN;                           %--Create U array with gp points and VFs
Uaver=zeros(gp,1);                               %--Create U vector to average the step velocities
%%------------------
ypperturbated=perturb(yplus,Np);                 %--perturbation function in linear space
ypro=round(ypperturbated);                       %--Round pertrubated positions
[sorty,sortu]=sortyv(ypro,uplus);                %---function to sort y+ and u+ vectors increasing monotically
dyvf=abs(diff(sorty))>=1;                         %--make sure that distance between vf is >1

   while (dyvf==false & Np> 0)                     %--Don't allow repeated elements
         ypperturbated=perturb(yplus,Np);
         ypro=round(ypperturbated);
         [sorty,sortu]=sortyv(ypro,uplus);
         dyvf=abs(diff(sorty))>=1;
   end
%%%----------Fill the step profile
lp=0;                                               %--Number of points to the left of the centroid
  for j=1:nypl-1
      %----Storage index of current and next step
      lp=floor(dy/2);
      if (sorty(j)-lp)<=1
          ini1step=1;
      else
         ini1step=sorty(j)-lp;                      %----Initial index 1st step  (Current step)
      end
      fni1step=ini1step+(dy);                       %----Final index 1st step
      ini2step=ini1step+(dy);                       %----initial index 2nd step
      if (sorty(j+1)-lp)<=1
          fni2step=1;
      else
          fni2step=sorty(j+1)-lp;                   %----Final index 1st step
      end

      if j==1                                       %----Fill the first vortical fissure
         U(ini1step:ini1step+lp,j)=linspace(ul1vf,sortu(j),lp+1);                                      %---Fill left edge
         U(ini1step+lp:fni1step,j)=linspace(sortu(j),0.5*(sortu(j)+sortu(j+1)),lp+1);                  %---Fill right edge
      else                                          %----Fill all the VF except the first one
         U(ini1step:ini1step+lp,j)=linspace(0.5*(sortu(j-1)+sortu(j)),sortu(j),lp+1);                  %---Fill left edge
         U(ini1step+lp:fni1step,j)=linspace(sortu(j),0.5*(sortu(j)+sortu(j+1)),lp+1);                  %---Fill right edge
      end
      %----Check for overlapping
      if ini2step > fni2step                        %---overlapping
         countovf=countovf+1;
      else
         U(ini2step:fni2step,j)=0.5*(sortu(j)+sortu(j+1));              %---not overlapping
      end   
  end   

  %-----Initial step and non slip boundary condition
  if (sorty(1)-lp)>1
     U(1:2,1)=linspace(0.5,ul1vf,2);
     U(2:sorty(1)-lp,1)=ul1vf;
%     U(1:sorty[1]-lp,1)=zeros(sorty[1]-lp,1);                           %---Zero step at the wall
  end
U=[zeros(1,nypl);U];                                %---Add not slip boundary condition
%%--------------------------------------------------
Uaver=mean(U,2,'omitnan');                           %---Average step velocity profiles

  %-----Function to sort vortical fissure y+
  function [ysort,usort]=sortyv(yraw,uraw)
     %----Input
     %----yraw   raw vector y without sort  (pertubated yp)
     %----uraw   raw vector u without sort
     %----Output
     %----ysort  sorted vector y in ascending order
     %----usort  sorted vector u according with original indexes in yraw
     [ysort,index]=sort(yraw,'ascend');
     usort=uraw(index);
  end

end
