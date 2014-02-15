!**********************************************************************!
ENTRY QESNuc_Get_TgtName(iTarget,TgtName)
!----------------------------------------------------------------------!
    TgtName=NucName(iTarget)
    QESNuc_Get_TgtName=.true.
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtAWght(iTarget)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtAWght=NucParam(2,iTarget)
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtMass(iTarget)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtMass=NucMass(0,iTarget)
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtNumb(TgtName)
!----------------------------------------------------------------------!
    QESNuc_Get_TgtNumb=-999
    do Target=-1,18
        if(NucName(Target)==TgtName)then
            QESNuc_Get_TgtNumb=Target
            return                                                     !because of this |
        endif                                                          !                |
    enddo                                                              !                V
    if(QESNuc_Get_TgtNumb==-999)then                                   !we never will reach this line with QESNuc_Get_TgtNumb.ne.-999, isn't it?!
        write(*,*) 'QESNuc_Get_TgtNumb ERROR: Unknown element!'        !why not stopping just here? also -999 is needed only for this
    endif
    return

!**********************************************************************!
ENTRY QESNuc_Get_TgtNcln(iNuAnu,iTarget)
!----------------------------------------------------------------------!
    Z=NucParam(1,iTarget)
    A=NucParam(2,iTarget)
    if(iNuAnu==1)then
        QESNuc_Get_TgtNcln=A-Z
    else
        QESNuc_Get_TgtNcln=Z
    endif
    return

!**********************************************************************!
ENTRY QESNuc_Print_TgtParam(iTarget)
!----------------------------------------------------------------------!
    write(*,'(A40)') '****************************************'
    write(*,'(A1,1X,A9,I2,1X,A1,2X,2A2,2(A3,1X,I3,1X),A1,1X,A1)') &
    '*','Nucleus #',iTarget,': ',adjustr(NucName(iTarget)),' (',&
    'Z =',int(NucParam(1,iTarget)),'A =',int(NucParam(2,iTarget)),')','*'
    write(*,'(A1,11X,A27,A1)') '*','---------------------------','*'
    write(*,'(A1,10X,A1,6X,A2,5X,A1,5X,A3,5X,A1)') '*','|','Nu','|','Anu','*'
    write(*,'(A1,1X,A36,1X,A1)') '*','------------------------------------','*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.4),1X,A1)') '*','E_nuBin','|',NucParam(3,iTarget),'|',NucParam(4,iTarget),'*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.4),1X,A1)') '*','P_Fermi','|',NucParam(5,iTarget),'|',NucParam(6,iTarget),'*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.4),1X,A1)') '*','T_Fermi','|',NucParam(7,iTarget),'|',NucParam(8,iTarget),'*'
    write(*,'(A1,1X,A36,1X,A1)') '*','-------------NucVolFact-------------','*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.8),1X,A1)') '*','e','|',NucVolFact(1,1,iTarget),'|',NucVolFact(1,2,iTarget),'*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.8),1X,A1)') '*','mu','|',NucVolFact(2,1,iTarget),'|',NucVolFact(2,2,iTarget),'*'
    write(*,'(A1,2X,A7,2(1X,A1,1X,F11.8),1X,A1)') '*','tau','|',NucVolFact(3,1,iTarget),'|',NucVolFact(3,2,iTarget),'*'
    write(*,'(A40)') '****************************************'
    QESNuc_Print_TgtParam=.true.
    return
