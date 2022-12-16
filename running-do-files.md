# How to make the .do files run
1. Clone the repository

2. Download the data from [here](https://www.openicpsr.org/openicpsr/project/113005/version/V1/view) (you can connect with your google account)

3. Complete the cloned repository with the data and log files downloaded in step 2:
   - 5 "Log-Files" folders, that go in `DataClean`, `DescStats`, `Figures`, and `Models`
   - The whole `Datasets` folder goes in `Work`
   - `Bootstrap-Files` and `Estimates` go in `Models`

4. Start x2goclient. In the parameters of your PSE session, give x2go access to your local repo (on your computer), following [these instructions](https://wiki.x2go.org/doku.php/doc:howto:x2goclient-file-sharing).

5. Start the session. Open the file explorer. Go to `~/media/disk`. There, you should find a folder with a long name that ends in "robust-micrometrics". This is the local repo!

6. Open the local repo, and go to `original-project/voxA`. Create a copy of the `voxAdir.ado` file and put it in a new folder called `ado`, that you have to create at the root of your session (`/home/j.doe`).

7. Open this `/home/j.doe/voxAdir.ado` file with a basic text editor (right-click > Éditeur de texte)

8. Line 5 should be something like this:

   ```stata
   	cd "/home/r.turquier/media/disk/_Users_remi_Documents_GitHub_robust-micrometrics/original-project/voxA/Work/`ddddd'"
   ```

   Update this path such that it points to the `Work` folder in your local repo (you'll need to change `r.turquier` and `_Users_remi_Documents_GitHub_`).

9. Change the paths at the beginning of each of the 4 "include" do files (located in the "Models", "Figures", "DescStats" and "DataClean" folders).

10. Open Stata and try running the do files you want to run. If it doesn't work and you don't know why, curse about how Stata is badly designed and not adapted to making reproductible open science. Then, ask for help :)
