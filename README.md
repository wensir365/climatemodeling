Climate Modeling
----------------

[ Introduction ]

This project, hosted by GitHub, is aimed to provide some demo code to those students who chose the course "Climate Modeling" for better understanding the necessary concepts and skills of numerical modeling. All the code are subject to GNU/GPL license and open to the public.

We would like to thank all the students of our department.


[ HowTo ]

Set Up Git on Local

- git config --global user.name "wensir365"					// your account ID in github.com
- git config --global user.email "wensir365@163.com"				// your email
- git config --global credential.helper cache					// turn on password cache
- git config --global credential.helper 'cache --timeout=3600'			// set time (in sec) for password cache
- git config --global http.sslVerify false					// in case you met "SSL" error

Clone and Commit

- git clone https://github.com/wensir365/climatemodeling.git			// get a copy of current repo
- cd climatemodeling								// way in
- ls -a										// want to take a look?
- cp ../climatemodeling-code/01.lorenz/ . -rp					// creat or copy something in
- git add 01.lorenz/								// "add" files in local repo
- git status									// want to take a look on local repo?
- git commit -m 'first try'							// "commit" in local repo with a message
- git remote add origin https://github.com/wensir365/climatemodeling.git	// set remote "add" point
- git push origin master							// push your commit

Xinyu Wen
Aug-3-2012
