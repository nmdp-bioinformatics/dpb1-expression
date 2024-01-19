import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { HlaMatchComponent } from './hla-match.component';

describe('HlaMatchComponent', () => {
  let component: HlaMatchComponent;
  let fixture: ComponentFixture<HlaMatchComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ HlaMatchComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(HlaMatchComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
