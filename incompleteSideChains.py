import chimera
import Midas

from chimera import replyobj
from Rotamers import useBestRotamers
from DockPrep.prefs import prefs, defaults, INCOMPLETE_SC

mols = chimera.openModels.list(modelTypes=[chimera.Molecule])
incompleteSideChains = "rotamers"
rotamerPreserve = True
rotamerLib=defaults[INCOMPLETE_SC]
replyobj.status("Mutating incomplete side chains")
targets = []
for m in mols:
        for r in m.residues:
                tmplRes = chimera.restmplFindResidue(r.type,
                                                False, False)
                if not tmplRes:
                        continue
                t_amap = tmplRes.atomsMap
                if 'CA' not in t_amap or 'CB' not in t_amap:
                        continue
                r_amap = r.atomsMap
                if 'CA' not in r_amap:
                        continue
                todo = ['CB']
                seen = set(['N', 'CA', 'C'])
                incomplete = False
                while todo:
                        aname = todo.pop()
                        seen.add(aname)
                        if aname not in r_amap:
                                incomplete = True
                                break
                        ta = t_amap[aname]
                        for n in ta.bondsMap.keys():
                                if n.name in seen:
                                        continue
                                if n.element.number == 1:
                                        continue
                                todo.append(n.name)
                if not incomplete:
                        continue
                if incompleteSideChains == "rotamers":
                        targets.append(r)
                        continue
                # mutate to gly/ala
                if 'CB' not in r_amap:
                        replyobj.info("Mutating %s (incomplete "
                                "side chain) to GLY\n" % str(r))
                        r.type = 'GLY'
                        continue
                replyobj.info("Mutating %s (incomplete side"
                                " chain) to ALA\n" % str(r))
                r.type = 'ALA'
                seen = set()
                for bbName in ['N', 'CA', 'C']:
                        try:
                                alist = r_amap[bbName]
                        except KeyError:
                                continue
                        for a in alist:
                                seen.add(a)
                todo = []
                for cb in r_amap['CB']:
                        todo.append(cb)
                deathRow = set()
                while todo:
                        a = todo.pop()
                        seen.add(a)
                        for nb in a.neighbors:
                                if nb in seen:
                                        continue
                                if nb.residue != a.residue:
                                        continue
                                todo.append(nb)
                                deathRow.add(nb)
                for dr in deathRow:
                        dr.molecule.deleteAtom(dr)
if incompleteSideChains == "rotamers":
        if targets:
                replyobj.info("Residues with incomplete"
                                        " side chains:\n")
                for t in targets:
                        replyobj.info("\t" + str(t) + "\n")
                replyobj.info("Replacing each by 'swapaa"
                        " same (residue atom spec)")
                replyobj.info(" lib %s" % rotamerLib)
                replyobj.info(" preserve %s'\n"
                                        % rotamerPreserve)
                useBestRotamers("same", targets, lib=rotamerLib,
                        ignoreOtherModels=False, preserve=rotamerPreserve)
        else:
                replyobj.info("No incomplete side chains\n")

Midas.write(mols, None, './temp.pdb')