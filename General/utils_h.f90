module utils_h

!Till: included whichn function
!2008-08-12

!Till: inluded wrapper for pause command
!2008-08-08



contains
INTEGER FUNCTION GetNumberOfSubstrings(string)
!@ This function returns the number of substrings within a string. @
!@ Valid substring separators are one or more blanks " " or tabs  @
!@ characters (char(9)).                                          @

implicit none
!INTENT(IN)
character(len=*),intent(in):: string
!LOCAL
integer:: i,n
!CODE
!If string length is zero
if (len(string) .eq. 0) then
  n= 0
!If string has only one character
else if (len(string) .eq. 1) then
  if ((string .eq. " ").or.(string .eq. char(9))) then
    n= 0
  else
    n= 1
  end if
!If string has more than one character
else
  n= 0
  !Account for valid character at first position
  if ((string(1:1) .ne. " ") .and. (string(1:1) .ne. char(9))) then
    n= n + 1
  end if
  !Account for subsequent substrings
  do i=1,(len(string)-1)
    if (((string(i:i) .eq. " ") .or. (string(i:i) .eq. char(9))) .and. &
       ((string(i+1:i+1) .ne. " ") .and. (string(i+1:i+1) .ne. char(9)))) then
      n= n + 1
    end if
  end do
end if
GetNumberOfSubstrings= n
END FUNCTION GetNumberOfSubstrings



INTEGER FUNCTION id_ext2int(ext_id, id_array)

	! converts external IDs (used in the input files) to internal numbering scheme
	! 2005-06-29, Till Francke
	!ii pointer statt array übergeben (loc)


	IMPLICIT NONE

	INTEGER, INTENT(IN):: ext_id		!external ID that is to be converted into internal numbering scheme
	INTEGER, DIMENSION(:), INTENT(IN):: id_array		!lookup array containing the external IDs in the order corresponding to their internal numbering scheme

	INTEGER	:: i		!auxiliary variables


	DO i=1,size(id_array)					!search entire array for occurence of ext_id
		if (id_array(i)==ext_id) then
			id_ext2int=i					!external ID found, return internal position
			return
		end if
	END DO

	id_ext2int=-1					!no matching external ID found
	return

END FUNCTION id_ext2int


function which1(boolarray, nowarning) result(indexarray1)
!returns index of the FIRST element that is TRUE in the input array
!if there are more, issue warning, if enabled
implicit none
logical,dimension(:),intent(in):: boolarray
logical, optional :: nowarning
integer,dimension(max(1,count(boolarray))):: indexarray
integer:: indexarray1
integer:: i
	if (count(boolarray) .eq. 0) then
	  indexarray(:)= 0
	else
	  indexarray(:)= pack((/(i,i=1,size(boolarray))/),boolarray)
	end if
	
	if (size(indexarray)>1 .AND.  (.NOT. present(nowarning))  ) then 
		write(*,*)'Ambiguous result in which1 (replicates of soil/veg combination in TC?). Saved init_conds questionable.'
		indexarray1=indexarray(1)		!use first occurence anyway
		!stop 
	else
		indexarray1=indexarray(1)
	end if
end function which1

function whichn(boolarray, return0) result(indexarray)
!returns indices of the ALL elements that are TRUE in the input array
!return0: determine, if empty array or 0 is returned when entire array is false
implicit none
logical,dimension(:),intent(in):: boolarray
integer,intent(in):: return0
integer,dimension(max(count(boolarray), return0)) :: indexarray
integer:: i

	if (count(boolarray) == 0) then
		if (return0 == 1) then
			indexarray= 0	 !return 0 when entire boolarray is false
		end if
	else
	  indexarray(:)= pack((/(i,i=1,size(boolarray))/),boolarray)
	end if
end function whichn



subroutine pause1
!wrapper for pause command (doesn't work on all platforms, so disable here)
	!pause
	return
end subroutine


function locase(string) result(lower)	!convert to lower case
character(len=*), intent(in) :: string
character(len=len(string)) :: lower
integer :: j
do j = 1,len(string)
  if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
       lower(j:j) = achar(iachar(string(j:j)) + 32)
  else
       lower(j:j) = string(j:j)
  end if
end do
end function locase




END MODULE utils_h