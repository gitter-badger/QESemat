!**********************************************************************!
FUNCTION MA_QES_Init(iCorV,idMA,MA_QES)
!----------------------------------------------------------------------!
!MA_QES(_eff) setter
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!

implicit none

logical:: &
    MA_QES_Init
real:: &
    MA_QES_eff

real,parameter:: &
    MA0 =1.015,   E0 =0.331,&
    MAl =MA0-0.026,El =E0-0.056,&
    MAh =MA0+0.026,Eh =E0+0.06,&
    MAll=MA0-0.031,Ell=E0-0.066,&
    MAhh=MA0+0.032,Ehh=E0+0.072,&
    MAc0=1.012,&
    MAcl=MAc0-0.031,MAcll=MAc0-0.06,&
    MAch=MAc0+0.031,MAchh=MAc0+0.061
integer &
    iCorV,idMA,Target
real &
    E_nu,MA_QES,&
    MA(-2:2)/MAll,MAl,MA0,MAh,MAhh/,E(-2:2)/Ell,El,E0,Eh,Ehh/,&
    MAc(-2:2)/MAcll,MAcl,MAc0,MAch,MAchh/

save
integer &
    CorV,dMA
real &
    cff_MA,MA_QES_const/0./

!cff_MA setting--------------------------------------------------------!
    CorV=iCorV
    dMA =idMA
    cff_MA=1.+dMA*.03
!----------------------------------------------------------------------!
    if(CorV==1)MA_QES_const=cff_MA*MA_QES
    if(CorV==3)MA_QES_const=MAc(dMA)
    MA_QES_Init=.true.
    return

!**********************************************************************!
ENTRY MA_QES_eff(Target,E_nu)
!----------------------------------------------------------------------!
    MA_QES_eff=MA(dMA)
    if(CorV==2)then
        if((Target.ne.0).and.(Target.ne.-1))then
            MA_QES_eff=MA(dMA)*(1.+E(dMA)/E_nu)
        endif
    else
        MA_QES_eff=MA_QES_const
    endif
    return
!----------------------------------------------------------------------!
endFUNCTION MA_QES_Init
!**********************************************************************!
