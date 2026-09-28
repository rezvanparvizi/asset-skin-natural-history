# Server notes — mininubio, screen, Positron, git

Practical how-to for this specific setup. Written to be followed
literally the first time and skimmed thereafter.

- [The setup, end to end](#the-setup-end-to-end)
- [One-time: SSH key for the server](#one-time-ssh-key-for-the-server)
- [One-time: SSH key for GitHub, from the server](#one-time-ssh-key-for-github-from-the-server)
- [One-time: clone and wire up the repo](#one-time-clone-and-wire-up-the-repo)
- [Positron: connecting to the server](#positron-connecting-to-the-server)
- [screen: the part nobody explains](#screen-the-part-nobody-explains)
- [Daily loop](#daily-loop)
- [Things that will bite you](#things-that-will-bite-you)

---

## The setup, end to end

```
U-M laptop (Positron)
      |
      |  VPN  (must be connected first, always)
      v
mininubio.ddns.med.umich.edu      <- ssh; no scheduler, jobs run in screen
      |
      |  the repo lives HERE, cloned from GitHub
      |  the data lives HERE, on large storage, never in git
      v
GitHub (private, personal account)   <- code only
```

You edit on the laptop through Positron's remote connection, but the
files, the R process, and the data are all on the server. Nothing heavy
crosses the VPN except your keystrokes and rendered plots.

---

## One-time: SSH key for the server

On the **laptop**, in a terminal (Positron has one: Terminal → New
Terminal):

```bash
ssh-keygen -t ed25519 -C "laptop-to-mininubio"
# accept the default path; set a passphrase if you want one
```

Copy the public key to the server:

```bash
ssh-copy-id <uniqname>@mininubio.ddns.med.umich.edu
# if ssh-copy-id is unavailable:
#   cat ~/.ssh/id_ed25519.pub | ssh <uniqname>@mininubio.ddns.med.umich.edu \
#     "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys"
```

Add a host alias so you can type `ssh mininubio`. On the **laptop**,
edit `~/.ssh/config`:

```
Host mininubio
    HostName mininubio.ddns.med.umich.edu
    User <uniqname>
    ServerAliveInterval 60
    ServerAliveCountMax 10
```

`ServerAliveInterval` keeps the connection from being dropped while
idle. It does **not** protect a running job — that is what `screen` is
for.

Test: `ssh mininubio`

---

## One-time: SSH key for GitHub, from the server

The server needs its own key, because the server is what talks to
GitHub. GitHub no longer accepts password authentication.

On the **server**:

```bash
ssh-keygen -t ed25519 -C "mininubio-asset"
cat ~/.ssh/id_ed25519.pub
```

Copy that output. In a browser on your laptop: GitHub → Settings → SSH
and GPG keys → New SSH key → paste → Add.

Verify, on the server:

```bash
ssh -T git@github.com
# expect: Hi <username>! You've successfully authenticated...
```

---

## One-time: clone and wire up the repo

On the **server**:

```bash
# 1. where code lives (backed up, small)
mkdir -p ~/projects && cd ~/projects
git clone git@github.com:<username>/asset-skin-natural-history.git
cd asset-skin-natural-history

git config user.name  "Your Name"
git config user.email "you@umich.edu"
git config core.fileMode false     # avoids spurious diffs on shared filesystems

# 2. where data lives (large). ASK which paths are backed up and which
#    are purged on a schedule before you put anything irreplaceable here.
df -h
BIG=/path/to/large/storage/$USER/asset      # <-- edit
mkdir -p "$BIG"/{data,results,figures,logs}
mkdir -p "$BIG"/data/{raw,objects,bpcells,bulk,clinical,external,inherited}

ln -s "$BIG/data"    data
ln -s "$BIG/results" results
ln -s "$BIG/figures" figures
ln -s "$BIG/logs"    logs

# 3. tighten permissions on the clinical directory
chmod 700 "$BIG/data/clinical"

# 4. check paths resolve
Rscript R/paths_check.R

# 5. make the job scripts executable
chmod +x jobs/*.sh

# 6. restore packages
Rscript -e 'install.packages("renv"); renv::restore()'
```

---

## Positron: connecting to the server

Positron is VS Code-based, so it uses the same Remote-SSH mechanism.

1. **Connect the VPN first.** Nothing below works without it.
2. In Positron, open the Command Palette (`Cmd/Ctrl+Shift+P`).
3. Run **Remote-SSH: Connect to Host** and pick `mininubio` (it reads
   your `~/.ssh/config`, which is why the alias was worth setting up).
   If the command is absent, install the Remote-SSH extension from the
   Extensions pane.
4. Once connected, **File → Open Folder** →
   `~/projects/asset-skin-natural-history`.
5. Positron's R console now runs **on the server**. Check it:

```r
Sys.info()[["nodename"]]   # should be the server, not your laptop
.libPaths()                # should be server paths
```

6. Connect Positron's git pane: it picks up the cloned repo
   automatically. Stage, commit, and push from the Source Control pane,
   or from the terminal — same thing.

**Interactive versus batch.** Use the Positron console for exploration:
inspecting objects, trying a plot, checking a marker. Use
`jobs/run.sh` for anything that takes more than a few minutes. The
console dies when the VPN drops; a screen session does not.

---

## screen: the part nobody explains

`screen` gives your process a terminal that keeps existing after you
disconnect. The confusion usually comes from the "attach / detach"
dance — which you can mostly skip.

### The pattern that avoids the dance entirely

```bash
jobs/run.sh analysis/05_sc_integration/05_1_harmony.R
```

Under the hood this is `screen -dmS <name>`, where `-d -m` means
**start it already detached**. The job begins, your prompt comes
straight back, and you never have to detach manually. Then:

```bash
tail -f logs/20260923_142301_03_1_harmony_freeze01_reference.log
```

Ctrl-C stops watching the log. It does **not** stop the job.

### The commands, for when you do need them

| Command | What it does |
|---|---|
| `screen -ls` | list sessions: name, PID, Attached/Detached |
| `screen -r <name>` | attach to a session |
| `Ctrl-A` then `d` | detach from a session, leaving it running |
| `screen -S <name> -X quit` | kill a session and its job |
| `screen -S <name>` | start a new *interactive* session |
| `Ctrl-A` then `[` | scrollback mode; arrows/PgUp to scroll, `q` to exit |

`Ctrl-A` is screen's prefix: press and release Ctrl-A, *then* press the
next key. Do not hold them together.

### Reading `screen -ls`

```
There are screens on:
        48213.03_1_harmony_reference    (Detached)
        48877.05_2_sweep_placebo        (Attached)
```

`Detached` = running, nobody watching. That is the normal, good state.
`Attached` = someone (probably another of your terminals) has it open.
If you get "there is no screen to be resumed matching X", the job
finished or crashed — check the log.

### Stale sessions

Sessions outlive their jobs if the shell is still alive. Clean up:

```bash
screen -wipe            # remove dead sessions from the list
```

### Interactive R inside screen

Occasionally useful for a long exploratory session you want to survive a
VPN drop:

```bash
screen -S explore       # new interactive session
R                       # start R inside it
# ... work ...
# Ctrl-A then d         -> detach; R keeps running
# later:
screen -r explore       # exactly where you left it
```

This is the one case where you detach by hand.

### Don't do this

- **Don't run heavy R on the login shell without screen.** A VPN blip
  kills a six-hour integration with no log.
- **Don't run two jobs writing to the same run directory.**
  `jobs/run.sh` refuses to start a second session with the same name,
  which covers the common case.
- **Don't rely on `nohup`** here. It works, but you lose the ability to
  attach and inspect, and `screen -ls` becomes your only job list.

---

## Daily loop

```bash
# laptop: connect VPN, open Positron, Remote-SSH to mininubio

# server, in the repo
git pull                                   # if you also edit elsewhere

# ... write or edit a script in Positron ...

git add analysis/07_sc_subcluster_fibroblast/07_2_sweep_resolution.R
git commit -m "05_2: resolution sweep 0.1-0.8 on full-cohort fibroblasts"

# commit BEFORE a run you intend to keep, so git_sha describes the code
jobs/run.sh analysis/07_sc_subcluster_fibroblast/07_2_sweep_resolution.R \
  --freeze freeze01 --cohort reference --labelset labelset01

jobs/status.sh                             # check on it
tail -f logs/<newest>.log                  # watch it

# when it finishes: look at the result, then record the decision
#   docs/decisions.md   <- why you chose what you chose
#   docs/runs.csv       <- written automatically by finalize_run()

git add docs/decisions.md
git commit -m "decisions: fibroblast resolution 0.3, rejected 0.5 (no markers)"
git push
```

Two habits worth building deliberately:

1. **Commit before a run you intend to keep.** `init_run()` warns on a
   dirty tree because a run whose `git_sha` does not describe its code
   is not reproducible.
2. **Write the `decisions.md` entry the same day.** Two minutes now
   versus an hour of reconstruction in six months.

---

## Things that will bite you

**VPN drop mid-command.** The `ssh` session dies; anything in `screen`
survives. Reconnect, `screen -ls`, carry on.

**Filling the disk with `.rds` files.** Every Layer 1 stage writes
multi-GB objects. Check `df -h` before big runs. Prefer BPCells on-disk
matrices over in-memory objects at this cell count.

**Purged scratch space.** If your large storage is a scratch filesystem,
it may be swept on a schedule. Confirm which paths are backed up before
your only copy of a long integration lives there.

**Memory, with no scheduler to protect you.** A 1.3M-cell object can
exhaust a shared machine and take other people's jobs down with it.
`jobs/status.sh` shows free memory; `jobs/run.sh` caps BLAS threads at
4 by default. Raise it only deliberately:

```bash
OMP_NUM_THREADS=16 jobs/run.sh <script>
```

**R version drift.** The server's R and Bioconductor will be upgraded at
some point. `renv.lock` is what makes freeze01's results reproducible
afterwards. Run `renv::snapshot()` at each freeze and commit the lock
file.

**BPCells absolute paths.** Objects copied from a colleague point at
*her* directories. See `docs/inherited_objects.md` and
`R/io.R::check_bpcells()`.

**`screen` vs `tmux`.** If the server has `tmux` and you prefer it, the
concepts map directly (`tmux new -d -s name`, `tmux ls`, `tmux attach
-t name`, `Ctrl-B d`). `jobs/run.sh` uses `screen` because that is what
this server is set up for.
