## Compiling and debugging sgpp on Windows using Visual Studio

(not VS CODE)

- Install [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- Install the following packages on the WSL: `sudo apt install g++ gdb make rsync zip`
- Fetch the remote headers as described [here](https://devblogs.microsoft.com/cppblog/intellisense-for-remote-linux-headers/)
- Edit the linux.props file 
- Start ssh on the WSL using `sudo service ssh start`
- Build and run