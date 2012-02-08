
module string_manip

! make everything in this module private
private

implicit none

! Define those subroutines that can be used publicly
public get_num_tokens, get_token

contains

!*******************************************************************************
!
! Subroutine: get_num_tokens
!
! Description: Determines how many whitespace-delimited words are in a string
!
!*******************************************************************************

subroutine get_num_tokens(string, token_num)

  implicit none

! Passed arguments

  character(*), intent(in) :: string

  integer, intent(out)     :: token_num

! Local variables

  integer :: string_loc  ! our location in the string
  integer :: iend        ! last non-whitespace character location

  string_loc = 1
  iend = len_trim(string)
  token_num = 0

  do while (string_loc .le. iend)

    if ( string(string_loc:string_loc) .le. ' ' ) then
      string_loc = string_loc + 1
    else

      do while ( string(string_loc:string_loc) .gt. ' ' )
        string_loc = string_loc + 1
      end do

      token_num = token_num + 1
    end if
  end do

end subroutine get_num_tokens

!*******************************************************************************
!
! Subroutine: get_token
!
! Description: Gets the num'th token in a string
!
!*******************************************************************************

subroutine get_token(string, num, token)

  implicit none

! Passed arguments
  
  character(*), intent(in)  :: string  ! The string to parse
  character(*), intent(out) :: token   ! The token to return

  integer, intent(in)       :: num     ! Which token to return

! Local variables

  integer   :: num_tokens
  integer   :: istart
  integer   :: iend
  integer   :: string_loc
  integer   :: token_count

  ! Uncomment the below chunk of code for a "safe" get_num_tokens at the 
  ! expense of calling get_num_tokens() each time a specific token is
  ! pulled from the string. When it's commented out, token will just be
  ! a blank string upon return

! call get_num_tokens(string, num_tokens)

! if (num .gt. num_tokens)
!   write(mdout, *) ' Error in get_token: Looking for more tokens than &
!                     &there are in string'
!   call mexit(6,1)
! end if

  ! Now get the num'th token

  token_count = 0
  istart = 1
  iend = len_trim(string)
  token = ' '

  do while (istart .le. iend)

    if (string(istart:istart) .le. ' ') then

      istart = istart + 1
      
    else

      do string_loc = istart, iend
        if ( string(string_loc:string_loc) .le. ' ' ) exit
      end do

      token_count = token_count + 1

      ! If this is the token we want, store it and return
      if ( token_count .eq. num ) then
        token = string(istart:string_loc-1)
        return
      end if

      istart = string_loc ! Move to the next token

    end if

  end do

end subroutine get_token

end module string_manip
