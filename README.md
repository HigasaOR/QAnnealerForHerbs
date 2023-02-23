# QAnnealerForHerbs


## File structure

* ```src/```: code
* ```data/```: where you put your data

## Prerequisites 

* For the following steps, if you see a @ sign at the end, it is not necessary to do this. But if you want to run the ```example()``` in an individual python script successfully, you'll have to.

1. Get the herb data package (`Structure_hunter_execution_pack.zip`). Unzip, and you will find everything in `/home/pikachu/www/CSCCP/Structure_hunter_execution_pack/`.
2. Place `scaffold_0511.structureIdx.bak` under `data/db/` @
3. Put the herbs scaffold data under `data/herbs`, i.e., copy everything in `scaffold_structure example/` in the herb data package to `data/herbs`. @
4. Put the herbs molecular weight data under `data/mw`, i.e., copy everything in `MW/` in the herb data package to `data/mw`. @

## How to Run

0. First go to the ```src/``` directory:
    ```
    cd src/
    ```

1. Run ```search_scaffold.py``` to find the corresponding ```.cIdx``` file of the target scaffold:
    ```
    python search_scaffold.py [sdf_file]
    ```
    You will see the searched result in terminal. Then, You should put the searched file under the ```data/cIdxs/``` directory. For example: if you see ```cIdx found: s0000000580``` in the terminal, you should put the ```s0000000580.cIdx``` file into ```data/cIdxs/```.

2. Generate DAU-acceptable problem inputs. For this part, visit ```dau_input.py``` to see how it's actually done.
    Please visit `sample_main.py` to see how are things combined. You can also run it:
    ```
    python sample_main.py [sdf_file] [cIdx_file] [mw_file] 
    ```
    Because we are still experimenting with DAU so you can design your own generating files.

3. Combine DAU settings and problem inputs. TBD.
