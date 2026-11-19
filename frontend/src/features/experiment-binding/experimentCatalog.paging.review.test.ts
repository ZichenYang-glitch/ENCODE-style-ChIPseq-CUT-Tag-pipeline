import { describe, expect, it, vi } from 'vitest';

const { list } = vi.hoisted(() => ({ list: vi.fn() }));
vi.mock('../../api/generated/catalog/catalog', () => ({ catalogListExperiments: list }));
import { loadSelectableExperiments } from './experimentCatalog';

const project = 'prj_' + '1'.repeat(32);
const sample = { sample_revision_id: 'smpr_' + '2'.repeat(32) };
const revision = (index: number) => 'expr_' + index.toString(16).padStart(32, '0');

describe('experiment revision pagination', () => {
  it('offers eligible experiments on every page through the final valid cursor', async () => {
    list.mockImplementation(async (_project: string, query: { cursor?: string; limit: number }) => {
      const offset = Number(query.cursor ?? 0);
      const size = Math.min(query.limit, 1001 - offset);
      return {
        project_id: project,
        experiments: Array.from({ length: size }, (_, i) => ({
          display_name: 'Experiment ' + (offset + i + 1),
          experiment_revision_id: revision(offset + i + 1),
          samples: [sample],
        })),
        next_cursor: offset + size < 1001 ? String(offset + size) : null,
      };
    });
    const choices = await loadSelectableExperiments(project);
    expect(choices).toHaveLength(1001);
    expect(choices[0].experimentRevisionId).toBe(revision(1));
    expect(choices.at(-1)?.experimentRevisionId).toBe(revision(1001));
    expect(list).toHaveBeenCalledTimes(21);
  });

  it('rejects a repeated experiment cursor', async () => {
    list.mockImplementation(async (_project: string, query: { cursor?: string }) => (
      query.cursor === 'again'
        ? { project_id: project, experiments: [], next_cursor: 'again' }
        : {
            project_id: project,
            experiments: [{
              display_name: 'Experiment 1',
              experiment_revision_id: revision(1),
              samples: [sample],
            }],
            next_cursor: 'again',
          }
    ));
    await expect(loadSelectableExperiments(project)).rejects.toThrow(
      'The experiment list repeated a page.',
    );
  });
});
