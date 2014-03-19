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
    MA0=9.64d-01,E0=3.56d-01,a0=0.846d+00,&
    a1=0.10583186025394,b1=0.36563826143285,c1=0.91005905942944,&
    a3=0.00245436835458,b3=0.45861634912994,c3=1.03324251687687
integer &
    iCorV,idMA,Target
real &
    E_nu,MA_QES,&
    a(-1:1)/a1,0.,a3/,b(-1:1)/b1,0.,b3/,c(-1:1)/c1,0.,c3/

save
integer &
    CorV,dMA
real &
    cff_MA,MA_QES_const/0./

!cff_MA setting--------------------------------------------------------!
    CorV=iCorV
    dMA =idMA
    cff_MA=1.+dMA*.045
!----------------------------------------------------------------------!
    if(CorV==1)MA_QES_const=cff_MA*MA_QES
    MA_QES_Init=.true.
    return

!**********************************************************************!
ENTRY MA_QES_eff(Target,E_nu)
!----------------------------------------------------------------------!
    MA_QES_eff=c(dMA)+b(dMA)/(E_nu+a(dMA))
    if(CorV==2)then
        if((Target==0).or.(Target==-1))then
            MA_QES_eff=MA0
        elseif(dMA==0)then
            MA_QES_eff=MA0*(1.+(E0/E_nu)**a0)
        endif
    else
        MA_QES_eff=MA_QES_const
    endif
    return
!----------------------------------------------------------------------!
endFUNCTION MA_QES_Init
!**********************************************************************!
