#define maximum(a,b)  (((a) > (b)) ? (a) : (b))

label nReinit = 10000;

dimensionedScalar alpha1Boundary;

volScalarField alpha1_ls = alpha1/sqrt(pow(alpha1,2)+pow(mag(fvc::grad(alpha1))*dimChange*deltaX,2));
surfaceScalarField faceSurfaceGradient(fvc::snGrad(alpha1));
scalarField fsg(mesh.nFaces(),0.0);

scalar G = 0.0;
scalar eps = 1.0e-5;
scalar del;
label it_stop = 0;
volScalarField alpha1_temp = alpha1;

for(int reinit=0; reinit<nReinit; reinit++)
{
    volScalarField alpha1_temp = alpha1;         
    forAll(alpha1,cellI)
    {                                  
        label fa[3];
        double an[3],bn[3],an2[3],bn2[3];

        fa[0]=prismFaceLabels[cellI][2];
        fa[1]=prismFaceLabels[cellI][3];
        fa[2]=prismFaceLabels[cellI][4];
          	      	      
	    for(int fI=0;fI<3;fI++)
        { 
                if(mesh.isInternalFace(fa[fI]))
                {
                    fsg[fa[fI]] = faceSurfaceGradient[fa[fI]];
                }
                else
                {
                    label iPatch=mesh.boundaryMesh().whichPatch(fa[fI]);
                    label iFace=mesh.boundaryMesh()[iPatch].whichFace(fa[fI]);
                    const surfaceScalarField& magSf2 = mesh.magSf();
                    if(magSf2.boundaryField()[iPatch].coupled())
                    {
                      fsg[fa[fI]] = faceSurfaceGradient[fa[fI]];                 
                    }
                    else
                    {
                      fsg[fa[fI]] = faceSurfaceGradient.boundaryField()[iPatch][iFace];
                    }                  
                }
                                       
                vector r = mesh.faceCentres()[fa[fI]]-mesh.C()[cellI];
              
                if((alpha1[cellI]>0.0 && r[0]>0.0) || (alpha1[cellI]<0.0 && r[0]<0.0))
                {  
                  an[fI]=min(0.0,fsg[fa[fI]]*mesh.Sf()[fa[fI]].x()/(mesh.magSf()[fa[fI]]+SMALL));
                  an2[fI]=pow(an[fI],2);
                }
                else if((alpha1[cellI]>0.0 && r[0]<0.0) || (alpha1[cellI]<0.0 && r[0]>0.0))
                {  
                  an[fI]=maximum(0.0,fsg[fa[fI]]*mesh.Sf()[fa[fI]].x()/(mesh.magSf()[fa[fI]]+SMALL));
                  an2[fI]=pow(an[fI],2);             
                }
                if((alpha1[cellI]>0.0 && r[1]>0.0) || (alpha1[cellI]<0.0 && r[1]<0.0))
                {  
                  bn[fI]=min(0.0,fsg[fa[fI]]*mesh.Sf()[fa[fI]].y()/(mesh.magSf()[fa[fI]]+SMALL));
                  bn2[fI]=pow(bn[fI],2);             
                }
                else if((alpha1[cellI]>0.0 && r[1]<0.0) || (alpha1[cellI]<0.0 && r[1]>0.0))
                {  
                  bn[fI]=maximum(0.0,fsg[fa[fI]]*mesh.Sf()[fa[fI]].y()/(mesh.magSf()[fa[fI]]+SMALL));
                  bn2[fI]=pow(bn[fI],2);              
                }
              
          }
          //Info << "-----------------------------------------------" << endl;
          scalar ai_max=an2[0];
          scalar bi_max=bn2[0];
          
          for (int i=0;i<3;i++) 
          {
             if(an2[i]>ai_max)    ai_max=an2[i];
             if(bn2[i]>bi_max)    bi_max=bn2[i];
          }
          G = Foam::sqrt(ai_max+bi_max)-1.0;
          
          alpha1_ls[cellI] = alpha1_ls[cellI]*G;
          //Info << "sign_ls[" << cellI << "] = " << sign_ls[cellI] << "  G = " << G << endl; 
                 
          alpha1[cellI] = alpha1[cellI] - deltaTau.value()*alpha1_ls[cellI];    
     }
    
     alpha1.correctBoundaryConditions();
     
     del = 0.0;
     forAll(alpha1,cellI)
     {
        del = maximum(del,mag(alpha1[cellI]-alpha1_temp[cellI]));
     } 
     if(del<eps)  break;
     it_stop = reinit;
}
Info << "Iteration = " << it_stop << "  del = "<< del << endl;
Info << "---------------------------------------" << endl;          
