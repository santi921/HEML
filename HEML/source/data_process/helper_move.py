import os

# crawl all sub directories of a folder up to depth of three and copy folder structure to new folder
root = "./compressed_frames/"
output = "./compressed_frames_bare/"

for root, dirs, files in os.walk(root):
    for dir in dirs:
        if not os.path.exists(os.path.join(output, dir)):
            os.makedirs(os.path.join(output, dir))
            # make a folder in each called embedding
            os.makedirs(os.path.join(output, dir, "embedding"))
            # make a folder in embedding called o
            os.makedirs(os.path.join(output, dir, "embedding", "o"))
            os.makedirs(os.path.join(output, dir, "embedding", "oh"))
            os.makedirs(os.path.join(output, dir, "embedding", "normal"))
            dir_in = os.path.join(root, dir)
            dir_out = os.path.join(output, dir)
            os.system(
                "cp {}/embedding/o/control {}/embedding/o/control".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/oh/control {}/embedding/oh/control".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/normal/control {}/embedding/normal/control".format(
                    dir_in, dir_out
                )
            )

            os.system(
                "cp {}/embedding/o/coord {}/embedding/o/coord".format(dir_in, dir_out)
            )
            os.system(
                "cp {}/embedding/oh/coord {}/embedding/oh/coord".format(dir_in, dir_out)
            )
            os.system(
                "cp {}/embedding/normal/coord {}/embedding/normal/coord".format(
                    dir_in, dir_out
                )
            )

            os.system(
                "cp {}/embedding/o/gradient {}/embedding/o/gradient".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/oh/gradient {}/embedding/oh/gradient".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/normal/gradient {}/embedding/normal/gradient".format(
                    dir_in, dir_out
                )
            )

            os.system(
                "cp {}/embedding/o/GEO_OPT_CONVERGED {}/embedding/o/GEO_OPT_CONVERGED".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/oh/GEO_OPT_CONVERGED {}/embedding/oh/GEO_OPT_CONVERGED".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/normal/GEO_OPT_CONVERGED {}/embedding/normal/GEO_OPT_CONVERGED".format(
                    dir_in, dir_out
                )
            )

            os.system(
                "cp {}/embedding/o/opt.xyz {}/embedding/o/opt.xyz".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/oh/opt.xyz {}/embedding/oh/opt.xyz".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/normal/opt.xyz {}/embedding/normal/opt.xyz".format(
                    dir_in, dir_out
                )
            )

            os.system(
                "cp {}/embedding/o/final.xyz {}/embedding/o/final.xyz".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/oh/final.xyz {}/embedding/oh/final.xyz".format(
                    dir_in, dir_out
                )
            )
            os.system(
                "cp {}/embedding/normal/final.xyz {}/embedding/normal/final.xyz".format(
                    dir_in, dir_out
                )
            )
