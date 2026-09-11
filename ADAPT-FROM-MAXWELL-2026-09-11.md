# Boltzmann — adapting Maxwell into the parallel library

*2026-09-11. Maxwell is the serial descendant of this code and has had a
year's worth of design work done to it. This note maps that work onto the MPI
pipeline: what ported unchanged, what changed shape because the library is
distributed, what does not port at all, and what the port found.*

*Every number below was measured in this session against MPICH 4.2.3, g++
13.3, at 1, 2, 3, 4, 5 and 7 ranks, with the old and new libraries built from
the same sources and driven by the same programs.*

---

## 0. The short answer

Maxwell's work divides cleanly in three.

**Ports unchanged, and is worth more here than there.** The graph
representation (edge arena, sorted-array adjacency), the honest memory
accounting, header hygiene, the `Jacobian` value type, `free_jacobian`, the
quiet `harvest`, `get_partitions`/`get_cost`, the no-tape guards, and the
build tree. All of it is done, measured, and verified. **19x fewer
allocations, 2.5x fewer bytes, 4.3x faster at one rank and 3.2x at two, with
a bit-identical Jacobian.**

**Ports, but means something different.** `RUN_TO_END`. In Maxwell the
argument is that a throw out of the middle of an operator strands whatever the
section allocated. Here the throw is also **not collective**, which is a
different and worse failure.

**Ports in two thirds.** Maxwell's `Tape` work is three things: the tape
*object* (done — §6.1), the library owning the checkpoint loop (done —
`run_tape`), and the public RAII surface, which assumes one object owns one
tape. Here a `Tape` on p ranks is p views of one distributed tape and must
make *collective* misuse unrepresentable, not just local misuse. That one is
designed in §6.1 and deliberately not written.

**Does not port.** Templating the scalar for second derivatives (§6.2), and
binomial checkpointing (§6.4), which interacts with the pipeline in a way
Maxwell's design note does not have to consider.

**What the port found, which is the real return on it.** Three defects that
exist only because the library is parallel, all of them in the *contract*
rather than the arithmetic, all of them previously undocumented:

1. **The checkpoint loop is not collective.** At 8 partitions over 3 ranks,
   rank 0 runs three passes of the caller's section and ranks 1 and 2 run
   four. Measured.
2. **The library does blocking MPI from inside the caller's operators**, so a
   section containing any collective deadlocks — confirmed under gdb with
   both stacks.
3. **A process could only ever open one tape**, because `initialize()` called
   `MPI_Init` and `finalize()` called `MPI_Finalize`. That is why this tree
   had no assertion suite: a test that checks 48 elementary functions needs
   48 tapes.

---

## 1. What the two libraries are

Both record a **linearised** graph rather than a value tape and compute the
Jacobian by **vertex elimination**. Both run the user's program more than
once: a profiling pass that counts vertices and edges without building them,
then productive passes that build one memory partition each, with
`check_memory()` throwing `BreakException` when the partition budget is spent
and `checkpoint()` restoring the independents and replaying.

The difference is what happens to a partition once it is built.

- **Maxwell** keeps every partition in one process. `check_memory()` merges
  the partition it just recorded into `intmed_map` and moves on.
- **This library** distributes them. `top_owner_idx` counts partitions down
  from the total, rank *r* owns partitions *r*, *r+p*, *r+2p*, …, and when a
  rank has recorded its partition it **eliminates locally, serialises the
  surviving edges as `cji_t` records, and sends them one hop around a ring**
  to rank *r−1*, which merges them into its own `intmed_map` and continues.
  The last rank in the pipeline assembles the Jacobian.

So `Process::send`, `Process::build`, `local_vertex_elimination_normal`,
`vector_run`/`normal_run`, `is_proc`, `is_final_rank` and the `MPI_EDGE`
datatype have no counterpart in Maxwell at all, and everything in Maxwell
that touches the graph or the budget has a counterpart here.

---

## 2. Ported unchanged, and measured

### 2.1 The graph representation — the biggest win, and bigger here

Maxwell's DESIGN-NOTES §5 measured a tape using **2.44x** the memory the
library thought it did: `Vertex` was 120 bytes of which 96 was two
`std::map`s, and every edge was one heap allocation plus two red-black tree
nodes.

Ported in full: `EdgeArena` (block-allocated edges with a free list),
`AdjArena` (size-classed slabs for adjacency blocks), and `Adjacency`, a
sorted array binary-searched, keeping `std::map`'s `first`/`second` names and
its insert semantics so the elimination code reads unchanged.

| | before | after |
|---|---|---|
| `sizeof(Vertex)` | 120 | **72** |
| bytes per edge, as counted by the budget | 24 | **56** |
| bytes per edge, actually | ~112 + 3 malloc headers | **56** |
| allocations, 20x30 lattice, 2 MB budget | 343 630 | **17 932** |
| bytes allocated, same | 20 450 009 | **8 152 681** |
| allocations, same tape, 60 kB budget | 345 360 | **19 075** |
| wall clock, 40x60 lattice, np=1 | 0.792 s | **0.183 s** |
| wall clock, same, np=2 | 0.218 s | **0.068 s** |
| wall clock, same, np=4 | 0.149 s | **0.067 s** |

Jacobian bit-identical in every case.

**Why the parallel library gets more out of this than Maxwell does.** Two
reasons, neither of which applies to a serial tape.

- Maxwell builds a graph and tears it down once. Here every rank builds,
  eliminates, serialises and **tears down a whole graph per partition**,
  *P/p* times. The allocator churn the arena removes is multiplied by the
  number of partitions.
- `Process::send` walks every surviving vertex's in-list to pack `cji_t`
  records for `MPI_Send`. Against `std::map` that was a red-black tree
  traversal per vertex; against a sorted array it is a contiguous scan. The
  arena pays twice — once in elimination, once in the pack loop that feeds
  the wire.

### 2.2 Honest memory accounting (Maxwell SVEGP-24)

`check_memory()` counted `sizeof(Edge)` per edge and ignored the two adjacency
entries that reference it. Now `EDGE_BYTES = sizeof(Edge) + 2*sizeof(AdjEntry)`
and `VERTEX_BYTES = sizeof(Vertex)`, and `get_heap()` reports what the arenas
actually hold so the ratio cannot drift unnoticed.

**In Maxwell an under-counted budget costs replay time. Here it also decides
the work distribution**: the partition count sets how many partitions each
rank owns, how long the pipeline is, and how big every message on the ring is.
Load balance across ranks is downstream of this number being right. The two
errors were partly cancelling — the vertex got 40% smaller as the edge got
more expensive — so the partition count at a fixed budget moved only from 10
to 9 on the benchmark above; the *meaning* of the budget is what changed.

The arenas are also told the budget before anything is recorded (SVEGP-31), so
a tight budget actually makes a rank small instead of leaving arena capacity
floating above it.

### 2.3 Header hygiene

`inc/Active.hpp` said `using namespace boltzmann::internals;` and `inc/API.hpp`
said `using namespace boltzmann;`, both at global scope, and `API.hpp` pulled in
`Vertex.hpp`, `Edge.hpp`, `Process.hpp`, `<mpi.h>` and `<sys/time.h>`.

So `#include "boltzmann.hpp"` put the entire graph representation and the MPI
headers into the caller's translation unit. That is worse here than in
Maxwell for a specific reason: **`<mpi.h>` is not the caller's business.** A
program that only differentiates does not need MPI's declarations, and a
program that does use MPI should include `<mpi.h>` itself rather than acquire
it by accident from an AD library and then discover that the library's idea of
which MPI it was built against has to match.

`boltzmann.hpp` is now `Typedefs`, `Active`, `Jacobian`, `API`, `BreakException`
and nothing else. Three example `main.cpp` needed `using namespace boltzmann;`
added; the rest already had it. `BreakException` moved out of
`boltzmann::internals` into `boltzmann` — the caller is required to catch it, so it
was public in everything but name, and only reachable because a public header
leaked the internal namespace.

### 2.4 The `Jacobian` value type — worth more here than in Maxwell

Maxwell's SVEGP-07 was about the leak: `harvest()` hands back a caller-owned
`double**` and there was no way to release it. Both `free_jacobian()` and the
`Jacobian` value type are ported.

The parallel argument is different and stronger. **Exactly one rank assembles
the Jacobian.** With the `double**` interface that fact was invisible:
`harvest()` simply did not touch `A` on the other ranks, so

```cpp
double ** A;              // uninitialised
harvest(m,n,A);
use(A[0][0]);             // fine on one rank, garbage pointer on p-1
```

is undefined behaviour on every rank but one — silently, because the
single-rank run you debug on is the one case that works. **Every example in
this tree is written that way.** A `Jacobian` cannot be in that state: it is
empty on precisely the ranks that did not assemble it, `is_harvesting_rank()`
answers the question before `harvest()` is called, and `data()` is one
contiguous row-major block so broadcasting it is well defined.

### 2.5 Small things

- `largeint` is `std::uint64_t` (SVEGP-08). **And it crosses the wire**: the
  `MPI_EDGE` datatype was built from `MPI_UNSIGNED_LONG`, which is four bytes
  on LLP64, so a heterogeneous job would have disagreed about the width of a
  vertex index. It is `MPI_UINT64_T` now, and the committed type is
  `MPI_Type_create_resized` to `sizeof(cji_t)` so that an array of them
  strides correctly the day a field is added and tail padding appears.
- Message tags are enumerators, not `#define`s in a header.
- `harvest(..., print_out)` (SVEGP-06), and the `intmed_map.size()` line that
  Maxwell's audit left open as a known defect — it shouted straight through
  the quiet path — is gone. On p ranks it was p times as loud.
- `checkpoint()` printed `number of partitions = N` to stdout unconditionally,
  from inside the loop, on every rank. `get_partitions()` is the way to ask.
- `~Process()` and `destroy_graph()` (SVEGP-27). The old teardown lived only
  in `finalize()`'s non-profiling branch **and only on the final rank**, so
  every other rank leaked its whole graph on every run.
- `TapeState` (SVEGP-26): the `Process` singleton and eleven file-statics are
  one object, reached through a `thread_local` pointer, and the communicator
  went with them. `run_tape()` owns the checkpoint loop's try/catch. Both
  §6.1.
- Guards on every entry point with no tape open, and `regress` calls all of
  them and asserts none crash — Maxwell's SVEGP-26 lesson, that deleting a
  global and letting the compiler find the uses has a blind spot: a function
  that read a file-static that was simply zero now dereferences a null
  pointer, and that compiles.
- The build tree, and the layout around it: `inc/`, `src/`, `examples/`,
  `build/` — Maxwell's, directory for directory. `test/` became `examples/`,
  because the directories in it are demonstrations of the library rather than
  tests of it; the three that *are* tests (`fdcheck`, `invariants`,
  `regress`) are new and sit alongside them, which is exactly how Maxwell
  arranges the same split. One `build/`, examples discovered from the
  filesystem, `common.mk`, `make test`, `make sanitize`, and `make ranks`
  (see §5). Every object depends on every header — the old rule named one or
  two each, so editing `Typedefs.hpp` rebuilt almost nothing, which in this
  library means one translation unit thinking a `Vertex` is 120 bytes while
  another thinks it is 72.

---

## 3. Ports, but means something different: `RUN_TO_END`

Maxwell's SVEGP-32 adds a second way for a pass to end: instead of throwing
`BreakException` once the target partition is recorded, return normally and
let the caller's section run to completion passively. Same partitions, same
Jacobian, same elimination cost, about twice the passive work.

Maxwell's argument is that the throw leaves the section at a point that
depends on the budget and appears nowhere in the source, so everything the
section allocated before it is stranded — once per pass, and the number of
passes is the partition count, so tightening the budget multiplies the leak.

That argument holds here unchanged, and this session reproduced it twice by
accident: the two assertion suites written for this port leaked 49 432 bytes
in 167 allocations and 2 880 bytes in 9 respectively, entirely in the test
code, entirely from `new active[]` inside a section that was abandoned
mid-operator. Both are fixed — one with a destructor, one with `RUN_TO_END` —
and the fix is commented in each file, because the leak is the lesson.

`examples/hybrdj1` is the sharper demonstration. Its checkpoint loop has **no
`try`/`catch`** — the section is a bare call to `fcn()`. Under
`BREAK_ON_TARGET` it dies with `terminate called after throwing an instance of
boltzmann::BreakException`, on every rank, at any budget that chunks. Nothing in
the tree had ever run it. It now sets `RUN_TO_END` and passes at 1, 2, 3 and 4
ranks, which is the point: the contract said "wrap the section in try/catch"
and said it nowhere, and `RUN_TO_END` removes the requirement rather than
restating it.

**The parallel argument is separate, and it is in §4.**

---

## 4. What the port found

### 4.1 The checkpoint loop is not collective

A rank leaves `while(checkpoint(x,y))` as soon as *it* has recorded its last
partition. With *P* partitions over *p* ranks the ranks get ceil/floor shares,
so they run different numbers of passes:

```
P=8 partitions, budget 500
  np=3   rank 0 : 3 passes, 2 partitions
         rank 1 : 4 passes, 3 partitions
         rank 2 : 4 passes, 3 partitions
  np=5   rank 0 : 2 passes    rank 2..4 : 3 passes
```

Invisible while the section is pure arithmetic, which is what every example in
this tree does. The moment it contains anything the ranks must agree on — a
collective, an I/O step, a shared RNG draw — they make different numbers of
calls to it.

`set_pass_mode(PASSES_COLLECTIVE)` makes every rank run `ceil(P/p)` passes.
Both *P* and *p* are known to every rank without communicating (the profiling
pass is identical on all of them), so this costs **no message**, does not move
the Jacobian by one bit, and costs at most one extra *passive* execution of
the section on at most *p−1* ranks — `is_proc()` is false throughout, so
nothing is allocated, eliminated or sent. Implemented, default off, and
`regress` asserts that it equalises the pass counts while leaving the
partition counts alone.

### 4.2 The library does blocking MPI from inside the caller's operators

`PASSES_COLLECTIVE` and `RUN_TO_END` together are still **not sufficient** to
make a collective in the section safe, and finding out why is the most useful
thing this port produced.

`check_memory()` runs at the end of every recorded operation. When the budget
is spent it calls `normal_run()`, which blocks in `MPI_Probe` waiting for the
next rank's edge list. So one rank can be inside the library's `MPI_Probe`
while another is inside the user's `MPI_Allreduce`. Both stacks, captured
under gdb, two ranks, `RUN_TO_END` and `PASSES_COLLECTIVE` both on:

```
rank 0   PMPI_Probe
         <- Process::local_vertex_elimination_normal   Process.cpp:461
         <- Process::normal_run                        Process.cpp:681
         <- Process::check_memory                      Process.cpp:706
         <- Process::unary_op                          Process.cpp:1014
         <- boltzmann::unary_op                            API.cpp:259
         <- boltzmann::active::operator=                   Active.cpp:64
         <- main

rank 1   PMPI_Allreduce
         <- main
```

Duplicating the communicator — which `initialize()` does — prevents message
*crosstalk*. It cannot prevent this: the problem is not which communicator the
messages are on, it is that one rank is blocked in the library while the other
is blocked in the user.

So the contract is: **the section between `checkpoint()` calls must contain no
MPI.** That was always true and was written down nowhere. It is now written
down in `inc/Typedefs.hpp` beside `pass_t`, with the stacks.

Lifting it is a roadmap item (§6.3), not a setting.

### 4.3 A process could only open one tape

`initialize()` called `MPI_Init` unconditionally and `finalize()` called
`MPI_Finalize`. Two consequences, neither documented:

- **A second `initialize()` in a process aborts the job**, because `MPI_Init`
  after `MPI_Finalize` is an error. This is why no example differentiates more
  than one function, and it is why the tree had no assertion suite — a test
  that checks 48 elementary functions needs 48 tapes. *The missing tests and
  that line of code are the same fact.*
- **A host application that uses MPI itself had its MPI torn down.** An
  optimiser that asks this library for a Jacobian and then wants to
  `MPI_Allreduce` the step finds MPI already finalized.

MPI is now initialised on demand only if nobody else has, finalised at process
exit via `atexit` and only if this library started it, and the per-tape
communicator is freed by `finalize()` rather than leaked (a program opening
tapes in a loop exhausted the implementation's communicator context ids).
Existing programs, which end immediately after `finalize()`, behave exactly as
before.

---

## 5. The verification harness, and the one invariant Maxwell does not have

Three suites, ported from Maxwell's `fdcheck` / `invariants` / `regress`:

- **`fdcheck`** — 48 elementary functions and operator forms against central
  differences, each on its own tape, each at a budget tight enough to force
  several partitions, asserted on every rank. It also asserts that **exactly
  one** rank claims the Jacobian.
- **`invariants`** — three tape shapes (deep-and-narrow lattice,
  wide-and-shallow, long scalar chain) at three budgets.
- **`regress`** — one test per defect: the six from the 2026-09-11 review, the
  behaviours ported from Maxwell, and a misuse section that calls every public
  entry point with no tape open and asserts none of them crash. 41 checks.

**`invariants` is shaped differently from its Maxwell ancestor, and the reason
is a finding.** Maxwell pins a bit-exact Jacobian hash. Here, the Jacobian
*cannot* be bit-exact across rank counts:

```
np=1 vs np=4, lattice kernel, budget giving 18 partitions
    max relative difference   8.795e-17      (one ulp)
    np=1 is bit-identical to the unchunked answer
    np=2, 3, 4, 5, 7 all agree with each other
```

That is exactly what Maxwell's DESIGN-NOTES §5 predicts: accumulation across
eliminations is ordered by *elimination* order, and splitting the graph over a
pipeline changes which vertex is eliminated where, and therefore the order in
which fill-in accumulates into a surviving edge. It is reassociation of a
floating-point sum, not a lost or double-counted edge.

So there are two assertions, and conflating them would make one of them a lie:

- **agreement**, checked at every rank count: the chunked Jacobian must match
  the unchunked one to within 8 ulp. Worst observed: 1 ulp.
- **pinned**, checked at np=1 only: bit-exact hash and exact elimination cost,
  the refactoring change-detector. Skipped above one rank rather than pinned
  per rank count, because pinning seven rows times four rank counts would be
  pinning the schedule, and the schedule is allowed to change.

`make ranks` runs all three at 1, 2, 3 and 4.  They live in `examples/`
beside the demonstrations, as Maxwell's do. That target is the whole point
of the harness existing in a parallel library: `make test` runs at one rank
count and cannot see a violation of the only invariant that is specifically
about being parallel.

**Status:** suite green at 1, 2, 3, 4, 5 and 7 ranks; whole tree builds with
zero warnings under `-Wall -Wextra`; clean under AddressSanitizer and
LeakSanitizer with **no library leaks at all**. Six example directories that
call an API removed years ago are frozen in `backup_20260911/` with a README
saying what they were and how to port one back, which is what Maxwell does
with its own pre-checkpoint examples.  `examples/` now holds only directories
that build and run, so `make list` needs no exceptions section.  The suite
still supports a `SKIP` file — it is the right answer for an example being
worked on rather than one retired — but nothing uses it.

---

## 6. What did not port, and what I would do next

### 6.1 The `Tape` RAII surface — two thirds done, one third designed

**Status: `TapeState` and the library-owned loop are implemented. The public
`Tape` is designed below and deliberately not written.**

Maxwell's `Tape` is a constructor and destructor in place of
`initialize`/`finalize`, a `run()` that owns the try/catch, and a `Jacobian`
return. Underneath it is `TapeState`, one object holding the `Process` and
everything that used to be a file-static, reached through a `thread_local`
pointer. It divides into three pieces with very different price tags.

#### Done: `TapeState` (SVEGP-26)

This library had exactly the split Maxwell's SVEGP-26 removed — a `Process`
singleton plus **eleven** file-statics in `API.cpp` (`indep_size`, `dep_size`,
the four 2-D dimensions, `memory_available`, `run_counter`, and the two
shadow-copy buffers) — and the 2-D stride bug fixed on 2026-09-11 lived in the
second one, in code that had no business existing apart from the tape it
described. All eleven are `TapeState` members now, reached through one
`thread_local` pointer. `Process::instance`, `get_proc_instance()` and
`del_proc_instance()` are gone and the constructor is public, because there is
no singleton left to enforce.

**No public API change and no example change**, which is what let the existing
net verify the whole thing: all seven pinned invariants came out bit-identical
at 1, 2, 3 and 4 ranks, `fdcheck` 48/48, `regress` 50/50, clean under ASan and
LSan.

**The parallel half, which Maxwell does not have: the communicator moved onto
the tape too.** `mpi_boltzmann_comm` was a file-global duplicated from
`MPI_COMM_WORLD`. With one tape per process that is indistinguishable from a
member; with two, the second `initialize()` overwrote the first tape's handle
and a later `finalize()` freed it out from under the first tape's ring
traffic. It is a `Process` member now and `~TapeState()` frees it — the same
class of bug as the shadow copies, one layer down, and it would have become a
live bug the moment anything opened two tapes.

`finalize()` is now three lines: finalise the graph, clear the current-tape
pointer, `delete`. Previously it put eleven statics back one at a time, and
forgetting one poisoned the next tape — which is precisely how `run_target`
not being reset came to be a bug in Maxwell.

**What the `thread_local` buys, stated more narrowly than Maxwell can.**
Maxwell says it removes the structural reason the library could never be
thread-safe, and that the arenas and the elimination have not been audited for
it. Both true here. But an MPI library has a *second* prerequisite: two
threads each driving a tape would each be doing point-to-point MPI on their
own communicator, and that requires the runtime to have been initialised at
`MPI_THREAD_MULTIPLE`. This library calls `MPI_Init`, which requests nothing,
and both MPICH and Open MPI then provide `MPI_THREAD_SINGLE`. So per-thread
tapes need `MPI_Init_thread` **as well as** this change, and the honest
statement is that SVEGP-26 removed one of two structural barriers. Two tapes
in one thread — which `regress` M-06 exercises — need neither.

#### Done: the library owns the checkpoint loop

`run_tape(x, y, body)`, four overloads matching the four `checkpoint`
overloads, in `inc/Run.hpp`. It puts the `try`/`catch` inside the library, so
the third of the three requirements on the caller's section — that it be
exception-safe — is retired rather than restated.

This is not a style preference. `examples/hybrdj1` had **no** try/catch and
died with an uncaught `BreakException` on every rank the first time anything
ran it; two of the three assertion suites written for the Maxwell port leaked
(49 432 bytes in 167 allocations, and 2 880 in 9) through scratch arrays
stranded by that same throw. Three authors out of three got it wrong and all
three were me. `hybrdj1` is now the pilot: it uses `run_tape` and contains no
`catch` anywhere, and it passes at 1, 2 and 4 ranks under the default
`BREAK_ON_TARGET` — verified with `RUN_TO_END` removed, so it is `run_tape`
doing the work and not the break mode.

`run_tape` returns the pass count, which was previously unobservable and is
the first thing you want when a budget behaves oddly; `get_partitions()` is
the other half, and `regress` asserts `passes == partitions + 1`.

They are free functions, not members, because of the next section.

#### Designed, not written: the public `Tape`

Maxwell's `Tape` makes misuse unrepresentable for one process. The parallel
analogue has to make **collective misuse** unrepresentable, and that is a
different exercise — so here is the design, and the reason each part of it is
the way it is.

`Tape t(n, m, budget)` on p ranks is not p independent tapes. It is p views of
one distributed tape, and they are coupled in ways a single-process RAII type
never has to think about:

| coupling | consequence if broken | what the type must do |
|---|---|---|
| all ranks must agree on `n`, `m` and the budget | the partition count differs per rank, so rank *r* waits on a ring message rank *r+1* will never send — a hang, not an error | constructor is collective and **checks**: one `MPI_Allreduce` over the three numbers, and it refuses loudly on disagreement rather than deadlocking twenty minutes later |
| all ranks must construct in the same order | two tapes' communicators get crossed | `MPI_Comm_dup` is already collective, so MPI itself enforces this — the type only has to not hide it |
| all ranks must destroy in the same order | `MPI_Comm_free` mismatch | destructor is collective; a `Tape` must not be destroyed inside a rank-dependent branch |
| the section must run the same number of times on every rank | see §4.1 | `run()` sets `PASSES_COLLECTIVE` itself, because a type whose whole purpose is to remove footguns should not leave that one loaded |

So the shape is:

```cpp
class Tape
{
public:
  explicit Tape( largeint n , largeint m , largeint budget ,
                 MPI_Comm comm = MPI_COMM_WORLD );   // COLLECTIVE
  ~Tape();                                            // COLLECTIVE

  active *  independents( const double * x0 );
  active ** independents( largeint rows , largeint cols , const double * x0 );
  void      dependent_shape( largeint rows , largeint cols );

  template<class Body> largeint run( active  & y , Body body );
  template<class Body> largeint run( active  * y , Body body );
  template<class Body> largeint run( active ** y , Body body );

  void dependent ( const active & y );
  void dependents( const active * y );
  void dependents( active ** y );

  Jacobian harvest( bool print_out = false );   // empty except on one rank
  bool     harvesting() const;

  largeint partitions() const;   // this rank
  largeint cost() const;         // this rank
  largeint memory() const;
  largeint heap() const;

  void set_elim_mode ( elim_t  );   // not yet meaningful here; see 6.5
  void set_break_mode( break_t );
  void set_pass_mode ( pass_t  );
private:
  internals::TapeState * st_;
  Tape( const Tape & );
  Tape & operator=( const Tape & );
};
```

**Four decisions worth writing down, because each could reasonably go the
other way.**

1. **The communicator is a constructor argument, defaulting to
   `MPI_COMM_WORLD`.** The free functions hard-code `MPI_COMM_WORLD`, which
   means the library cannot be used by a program that has already split its
   ranks into a grid — and a PDE solver that wants a Jacobian per subdomain is
   exactly the caller this library is for. Taking the communicator is the
   single largest capability the `Tape` adds over the free functions, and it
   is free: `TapeState` already duplicates whatever it is given.

2. **The constructor validates collectively and the destructor is
   collective.** An `Allreduce` of three integers at tape construction is
   nothing next to a tape, and it converts the worst failure mode this library
   has — a hang whose cause is a typo in one rank's budget — into a message on
   rank 0. The cost is that `Tape` cannot be constructed inside a
   rank-dependent branch, which is a constraint worth having explicit rather
   than latent.

3. **`run()` turns on `PASSES_COLLECTIVE`.** It is the only mode under which
   the loop returns the same value on every rank (§4.1), the cost is at most
   one extra passive pass on at most p−1 ranks, and the Jacobian does not
   move. Leaving it off by default in the free functions is right — it changes
   observable pass counts, and existing programs should not change under
   them — but a new RAII surface has no existing programs to protect.

4. **`harvest()` returns an empty `Jacobian` off the assembling rank rather
   than broadcasting.** Broadcasting would be friendlier and is one line for
   the caller to add; doing it silently would make `harvest()` collective,
   which means a caller who harvests inside `if(rank==0)` — the natural thing
   to write — hangs. Between a friendly API that can deadlock and a blunt one
   that cannot, take the blunt one and document `data()`.

**What is deliberately *not* in it.** No `Tape::scratch()` handing out work
arrays: Maxwell floats the idea, and with `RUN_TO_END` and `run_tape` the leak
it addresses is already gone twice over. No nesting and no second derivatives
— that is §6.2 and needs the scalar templated. No attempt to make MPI in the
section safe; that is §6.3 and lives in the library, not the surface.

**Sequencing.** This should land *after* §6.3, not before. §6.3 changes when
the library communicates, `run()` is where that change is observable, and
doing the surface first means doing `run()` twice.
### 6.2 Templating on the scalar type

Unchanged from Maxwell's analysis: second derivatives need
`ActiveT<Active>`, not just multiple tapes, and it is the largest change by
far. Nothing about being parallel makes it easier or harder.

### 6.3 The blocking exchange — the parallel-only roadmap item

§4.2 is the thing to fix, and the fix has the same shape as Maxwell's
`RUN_TO_END`: **do not do library work at an arbitrary point inside the
caller's section.** Concretely, `local_vertex_elimination_normal`'s
`MPI_Probe`/`MPI_Recv` becomes an `MPI_Irecv` posted once per pass and drained
at the `checkpoint()` boundary, with `send()` becoming `MPI_Isend` into a
retained buffer. The cost is one buffered edge list per rank; the gain is that
the library's communication happens only where every rank is guaranteed to be,
which is the one place the loop is already collective under
`PASSES_COLLECTIVE`. That, plus `RUN_TO_END`, is what would let the contract
change from "no MPI in the section" to "collectives in the section are fine".

### 6.4 Binomial checkpointing (Revolve)

Maxwell's REVOLVE-DESIGN note observes that the library already implements the
c = 1 case — one snapshot, the independents at step zero, replayed for every
partition, `P(P+1)/2` partition-executions — and that c > 1 is the same
mechanism with more slots and a schedule.

That observation holds here, and the arithmetic is more favourable: the
pipeline already divides the elimination across p ranks, so what Revolve
attacks — the quadratic replay cost — is the part that does *not* parallelise.
Today's cost is `P(P+1)/2` executions of the section **on every rank**, since
every rank replays the whole program; Revolve would cut that to nearly linear
in P for all of them at once. It is probably the single largest remaining
speedup, and it is orthogonal to everything in §6.1–6.3.

The complication Maxwell does not have: a snapshot must restore not only the
caller's state but the counters that make vertex indices a deterministic
function of the operation sequence, **and those counters must come out
identical on every rank**, or partition *p*'s output vertices will not meet
partition *p+1*'s input vertices in `intmed_map`. That is checkable — it is
exactly what `make ranks` tests — but it should be designed before it is
built.

### 6.5 Smaller, all cheap

- `Vertex::eliminate()` returns a Markowitz-style cost that now feeds
  `get_cost()`; it could also drive elimination *ordering*, which nothing
  does.
- `harvest()` allocates `A` only on the assembling rank and leaves the
  caller's pointer untouched elsewhere. The `Jacobian` form fixes this for new
  code; the `double**` form should probably allocate a zeroed matrix
  everywhere, or be retired.
- Forward vs reverse elimination (Maxwell SVEGP-23) is selectable there and
  not here. Reverse is what the pipeline does. Forward in a pipeline is a
  genuinely different data flow — edges would travel the other way round the
  ring — so this is a design question, not a port.
- `doc/` is the one part of Maxwell's layout with no counterpart here:
  Maxwell keeps a LaTeX paper in it.  There is nothing to put in one yet, and
  an empty directory is not worth committing.
- The three `checkpoint` overloads still have copy-pasted bodies, and the
  2-D one is still where the stride bug lived. One implementation over a
  range would have made it impossible.
- No MPI return code is checked anywhere.


---

## 7. The name: Boltzmann

Adopted as a clean break on 2026-09-11. `namespace boltzmann`,
`boltzmann.hpp`, `libboltzmann.a`, `BOLTZMANN_*`. Nothing named `demon`
survives in the tree.

**Why.** Maxwell's demon is *one* gatekeeper deciding, molecule by molecule,
what may pass under an information budget — which is exactly what
`check_memory()` and `BreakException` do in the serial library, and why that
one is called Maxwell. Boltzmann is the canonical next step: from a single
demon to a statistical **ensemble**, many identical copies of the same system
considered at once. That is literally this library. Every rank runs the same
program and they differ only in which partition of the tape each one records,
which is an ensemble in the precise sense rather than a metaphor. And
"Maxwell–Boltzmann" is a pairing nobody needs explained, so "Boltzmann is the
parallel Maxwell" carries its own documentation.

**What else was considered, and why not.**

- **Bennett.** Charles Bennett resolved the demon paradox through the cost of
  erasure, *and* Bennett (1973) is the origin of the logarithmic checkpointing
  that Revolve implements. Almost too apt — it names §6.4 rather than the
  library, so it is worth holding until binomial checkpointing actually lands.
- **Gibbs.** The Gibbs ensemble is if anything a more literal fit for "many
  replicas of one system", and it is shorter. It collides with Gibbs sampling
  and the Gibbs phenomenon, both of which a scientific-computing reader meets
  first.
- **Szilard.** The engine that first tied information to entropy, but it is
  the *single*-molecule story, so it reads as another sibling rather than a
  scale-up.
- **Carnot.** The ring pipeline is a cycle and Carnot is about maximum
  efficiency, which is pretty, but Carnot is about heat engines rather than
  information, so it drops the thread Maxwell started.

### The rename was done by a script, and the script had to be fixed twice

`rename.sh` does it: dry-run by default, bounded to this tree's own sources,
word-bounded so an identifier that merely contains the letters is untouched,
and it reports rather than rewrites lines where "demon" is prose about the
thought experiment. It is idempotent.

Both bugs it had are worth recording, because both are the kind that makes a
rename script worse than a careful `sed`.

**`\bdemon\b` does not reach an identifier component.** After the first
`--apply` the namespace, the header and the archive were renamed and five
things were not: `-ldemon` (the boundary fails after `l`, which is a word
character), `mpi_demon_comm`, `mpi_started_by_demon` and `demon_mpi_shutdown`
(underscore is a word character too). The build failed at the link step with
`cannot find -ldemon`, which is the *good* case. Component rules were added
requiring an adjacent `_` or the `-l` prefix — which still leave the word
"demonstration" in this note alone, verified.

**The detector and the rules have to agree about what counts as an
occurrence.** The "what would change" scan used `grep -w`, so once the
word-bounded rules had run it found nothing left, printed **"nothing to
do"** and exited 0 — reporting success on a tree that no longer linked. It now
scans with the same breadth the rules cover, and the script ends with a
completeness check that greps for any remaining `demon` and exits non-zero
with the list. A half-rename that still *compiles* is the outcome to be afraid
of; that check is what makes it impossible.

### Verified after the rename

- Library builds with zero warnings under `-Wall -Wextra`.
- `make ranks` green: `invariants` 7 rows / 0 failures, `fdcheck` 48 checks /
  0 failures, `regress` 41 checks / 0 failures, at 1, 2, 3 and 4 ranks.
- `make test` green: all fourteen examples.
- `make sanitize` clean under AddressSanitizer and LeakSanitizer.
- Pinned Jacobian hashes unchanged, so the rename moved no arithmetic.
